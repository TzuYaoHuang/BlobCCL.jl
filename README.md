# BlobCCL

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://TzuYaoHuang.github.io/BlobCCL.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://TzuYaoHuang.github.io/BlobCCL.jl/dev)
[![CI workflow status](https://github.com/TzuYaoHuang/BlobCCL.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/TzuYaoHuang/BlobCCL.jl/actions/workflows/ci.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/TzuYaoHuang/BlobCCL.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/TzuYaoHuang/BlobCCL.jl)

BlobCCL.jl identifies and characterizes discrete blobs (bubbles or droplets) in volume-of-fluid fields (in n-D array form) produced by multiphase solver such as [InterfaceAdvection.jl](https://github.com/TzuYaoHuang/InterfaceAdvection.jl). It implements the Informed Component Labeling (ICCL) algorithm of [Hendrickson, Weymouth & Yue (2020)](https://doi.org/10.1016/j.compfluid.2019.104373), a flood-fill connected-component labeling scheme that uses local interface-normal information to correctly split components that are connected only through a thin, under-resolved neck, avoiding the spurious merging that plain connectivity-based labeling produces near the grid scale. For each labeled blob it returns the volume-fraction-weighted volume and centroid, from which an equivalent radius can be computed.
