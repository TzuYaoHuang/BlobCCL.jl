using BlobCCL


function generate_diffuse_droplet!(f::AbstractArray{T,2}, center::Tuple{T,T}, radius::T, thickness::T; ε=T(0.01)) where T
    for I in CartesianIndices(f)
        dist = sqrt((I[1] - center[1])^2 + (I[2] - center[2])^2)
        
        # Tanh profile: α approaches 0 inside the droplet, 1 outside, intermediate at the interface
        alpha = 0.5 * (1.0 + tanh((dist - radius) / thickness))

        if alpha>1-ε alpha = 1 end
        if alpha<ε alpha = 0 end
        
        # Use minimum to allow droplets to coexist without overwriting each other with background (1.0)
        f[I] = min(f[I], alpha)
    end
end


# 1. Initialize 100x100 domain with background fluid (α = 1.0)
Ng = (100, 100)
grid = ones(Float64, Ng)

# 2. Inject droplets with intermediate volume fractions at the boundaries
# Droplet 1: Center (30, 30), Radius 10.0, Interface thickness 1.5
generate_diffuse_droplet!(grid, (30.0, 30.0), 10.0, 1.5)

# Droplet 2: Center (70, 65), Radius 15.0, Interface thickness 2.0
generate_diffuse_droplet!(grid, (70.0, 65.0), 15.0, 2.0)

generate_diffuse_droplet!(grid, (80.0, 15.0), 8.0, 2.0)

# 3. Execute the algorithm
# isblob evaluates α < 1.0 - eps, capturing the core and the diffuse interface.
# blobval evaluates 1.0 - α, correctly accumulating fraction-weighted properties.
labels, blobs = BlobCCL.LabelAnalyzeBlob(grid; blobtarget=0)

# 4. Output validation
println("Domain Size: ", Ng)
println("Identified $(length(blobs)) distinct droplets.\n")

for b in blobs
    # Analytical volume of a 2D circle is π*r^2. 
    # The fraction-weighted numerical volume should closely match this.
    if b.label == 1
        analytical_vol = pi * 10.0^2
        expected_center = [30.0, 30.0]
    else
        analytical_vol = pi * 15.0^2
        expected_center = [70.0, 65.0]
    end

    println("Droplet $(b.label):")
    println("  Computed Volume:   $(round(b.volume, digits=4))")
    println("  Analytical Volume: $(round(analytical_vol, digits=4))")
    println("  Computed Centroid: [$(round(b.centroid[1], digits=4)), $(round(b.centroid[2], digits=4))]")
    println("  Expected Centroid: $expected_center")
    println("  Error in Volume:   $(round(abs(b.volume - analytical_vol)/analytical_vol * 100, digits=4))%")
    println("---")
end
