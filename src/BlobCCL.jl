module BlobCCL

using StaticArrays

include("Blob.jl")
export BlobData, BlobRadius

include("ICCL.jl")
export LabelAnalyzeBlob

end