module BlobCCL

using StaticArrays

struct BlobData{T,D}
    label::Int
    volume::T
    centroid::SVector{D,T}
end

BlobRadius(a::BlobData{T,2}) where {T} = sqrt(a.volume/π) |> T
BlobRadius(a::BlobData{T,3}) where {T} = cbrt(3a.volume/4π) |> T

export BlobRadius

# Assuming these are available in your environment.
import WaterLily: loc, inside, δ

function LabelAnalyzeBlob(f::AbstractArray{T,D}; blobtarget=0) where {T,D}
    Ng = size(f)
    labels = zeros(Int, Ng)
    current_label = 0

    Blobs = Vector{BlobData{T,D}}()
    queue = Vector{CartesianIndex{D}}(undef, prod(Ng))

    isblob(α) = abs(α - blobtarget) < 1 - 100eps(T)
    blobval(α) = abs((1 - blobtarget) - α)
    
    domain = inside(f)

    @inbounds for I ∈ domain
        if isblob(f[I]) && labels[I] == 0
            current_label += 1
            labels[I] = current_label

            queue[1] = I
            head, tail = 1, 1

            vol = zero(T)
            centroid = zero(SVector{D,T})

            while head <= tail
                curr = queue[head]
                head += 1

                val = blobval(f[curr])
                vol += val
                centroid += loc(0, curr) * val

                # Dimensionally generic neighbor traversal
                for i in 1:D
                    for s in (-1, 1)
                        n = curr + s*δ(i, curr)
                        
                        n ∉ domain && continue
                        
                        if isblob(f[n]) && labels[n] == 0
                            labels[n] = current_label
                            tail += 1
                            queue[tail] = n
                        end
                    end
                end
            end
            push!(Blobs, BlobData(current_label, vol, centroid / vol))
        end
    end

    return labels, Blobs
end

export LabelAnalyzeBlob

end