# Code for the paper: visualizing the wavenumber content of a point pattern

[![DOI](https://zenodo.org/badge/679169174.svg)](https://zenodo.org/badge/latestdoi/679169174)

To run the code, follow these steps:

1. Clone this repository
2. Open Julia and navigate to the repository ([install Julia](https://julialang.org/downloads/) if necessary)
3. Activate the project and instantiate the packages to the correct version by running
   ```julia
   import Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```
4. Build RCall if not done previously
   ```julia
   Pkg.build("RCall")
   ```
5. Run the file

See the [Julia documentation](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project) for more details on environments.

## Files

- `Fig1_introexample.jl`: Figure 1, a Thomas process example.
- `Fig2_lgcp.jl`: Figure 2, a log-Gaussian Cox process example.
- `Fig3_lansingwoods.jl`: Figure 3, analysis of the Lansing Woods data.
- `Fig4_impulseresponse.jl`: Figure 4, the impulse responses for certain regions.

## Additional code

The folder `src` contains a minimal implementation of the isotropic spectral estimation method proposed by Rajala et al. (2023).

## Note on R

Part of the simulation and the Lansing Woods data require the `R` package `spatstat`.
Note that `R` and its packages are not version controlled by the Julia environment.
For convenience, we provide the file  `session_info.jl` which can be run from Julia to print the version of `R` and `spatstat` which is installed.
A summary of the package versions we used is given below (this is the relevant subset of the output from running `sessionInfo()`, not the complete output)

```R
R version 4.3.0 (2023-04-21)
Platform: x86_64-apple-darwin20 (64-bit)
Running under: macOS Ventura 13.3.1

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] spatstat_3.0-6         spatstat.linnet_3.1-1  spatstat.model_3.2-4  
[4] rpart_4.1.19           spatstat.explore_3.2-1 nlme_3.1-162          
[7] spatstat.random_3.1-5  spatstat.geom_3.2-1    spatstat.data_3.0-1   

loaded via a namespace (and not attached):
 [1] mgcv_1.8-42           Matrix_1.5-4          lattice_0.21-8       
 [4] splines_4.3.0         abind_1.4-5           polyclip_1.10-4      
 [7] deldir_1.0-9          goftest_1.2-3         spatstat.sparse_3.0-1
[10] grid_4.3.0            compiler_4.3.0        spatstat.utils_3.0-3 
[13] tensor_1.5
```

## General implementation

Note that the methodology of the paper is in the package [PointProcessFilters.jl](https://github.com/SDS-EPFL/PointProcessFilters.jl)

## References
- Grainger, J. P., Rajala, T. A., Murrell, D. J., & Olhede, S. C. (2023). Visualizing the Wavenumber Content of a Point Pattern. arXiv preprint arXiv:2306.04198.
- T. A. Rajala, S. C. Olhede, J. P. Grainger, and D. J. Murrell, "What is the Fourier transform of a spatial point process?" IEEE Transactions on Information Theory, 2023.
- Jake P. Grainger. (2023). "SDS-EPFL/PointProcessFilters.jl: v0.1.0 (v0.1.0)." Zenodo. https://doi.org/10.5281/zenodo.8252035
