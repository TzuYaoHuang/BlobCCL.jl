struct BlobData{T,D}
    label::Int
    volume::T
    centroid::SVector{D,T}
end

BlobRadius(a::BlobData{T,2}) where {T} = sqrt(a.volume/π) |> T
BlobRadius(a::BlobData{T,3}) where {T} = cbrt(3a.volume/4π) |> T