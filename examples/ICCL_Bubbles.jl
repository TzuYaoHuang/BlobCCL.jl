using BlobCCL, InterfaceAdvection

# 1. Initialize 100x100 domain with background fluid (α = 1.0)
Ng = (100, 100)
grid = zeros(Float64, Ng)
alpha = grid*0
nhat = zeros((Ng...,2))

r1 = 10
r2 = 20
δ=0.9
d = r1+r2+δ
cen1 = (15,15)
cen2 = @. cen1 + d*(cosd(30),cosd(60))

sdf(xx) = min(√sum(abs2,xx.-cen1)-r1, √sum(abs2,xx.-cen2)-r2)

applyVOF!(grid,alpha,nhat,sdf)

# 3. Execute the algorithm
# isblob evaluates α < 1.0 - eps, capturing the core and the diffuse interface.
# blobval evaluates 1.0 - α, correctly accumulating fraction-weighted properties.
labels, blobs = BlobCCL.LabelAnalyzeBlob(grid; blobtarget=1)

# 4. Output validation
println("Domain Size: ", Ng)
println("Identified $(length(blobs)) distinct droplets.\n")

for b in blobs
    println("Droplet $(b.label):")
    println("  Computed Volume:   $(round(b.volume, digits=4))")
    println("  Computed Centroid: [$(round(b.centroid[1], digits=4)), $(round(b.centroid[2], digits=4))]")
    println("---")
end