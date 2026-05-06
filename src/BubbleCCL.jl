module BlobCCL

using StaticArrays

struct BlobData{T,D}
    label::Int
    volume::T
    centroid::MVector{D,T}
end

import WaterLily: loc, inside
import InterfaceAdvection: δd
function LabelAnalyzeBlob(f::AbstractArray{T,D}; multitarget=0) where {T,D}
    Ng = size(f)
    labels = zeros(Int, Ng)
    current_label = 0

    Blobs = Vector{BlobData{T,D}}()
    queue = Vector{CartesianIndex{D}}(undef, prod(Ng))

    isblob(α) = abs(α-multitarget) < 1-100eps(T)
    blobval(α) = abs((1-multitarget) - α)

    @inbounds for I∈inside(f)
        if isblob(f[I]) && labels[I] == 0
            current_label += 1
            labels[I] = current_label

            queue[1] = I
            head, tail = 1, 1

            vol = zero(T)
            centroid = zero(MVector{D,T})

            while head <= tail
                curr = queue[head]
                head += 1

                val = blobval(f[curr])
                vol += val
                centroid .+= loc(0,curr).*val

                for d∈[-3,-2,-1,1,2,3]
                    n = curr + δd(d,curr)
                    !insidecell(n) && continue
                    if isblob(f[I]) && labels[I] == 0
                        labels[n] = current_label
                        tail += 1
                        queue[tail] = n
                    end
                end
            end
            push!(Blobs, BlobData(current_label, vol, centroid./vol))

        end
    end

    return labels, Blobs
end

end