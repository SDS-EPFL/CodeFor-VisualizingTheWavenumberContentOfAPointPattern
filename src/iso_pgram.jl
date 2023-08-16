module IsoPP
    using SpecialFunctions, LinearAlgebra, Meshes

    """
        periodogram_iso(x::PointSet, region::Box, t,  taper_a::Real = 25)

    Computes the isotropic spectral estimate of a given collection of points.
    """
    function periodogram_iso(x::PointSet, region::Box, t,  taper_a::Real = 25)
        @assert taper_a > 0
        sl = sides(region)
        Vol = measure(region)
        lambda = length(x)/Vol
        
        # Compute periodogram
        value = iso_sum(x, t, sl, taper_a)
        bias = ppgram_iso_bias(sl, t, taper_a)

        lambda2 = (lambda^2-lambda/Vol) # remove bias in lambda squared estimate
        sdf = lambda .+ 2value .- lambda2 .* bias

        return (t=t, sdf=sdf, lambda=lambda)
    end

    ## bias
    """
        ppgram_iso_bias(sl, t, taper_a, n = (150, 200))

    Approximates the bias in the isotropic periodogram
    """
    function ppgram_iso_bias(sl, t, taper_a, n = (150, 200))    
        ## Numerical integration. A bit slow, need to optimsize.
        angs = range(0, pi, length = n[1])[2:end] # symmetric, half is enough
        da = angs[2] - angs[1]
        rmax = sqrt(sum(abs2, sl)) # afterwhich setcov = 0
        rgrid = range(0, rmax, length = n[2])[2:end]
        dr = rgrid[2] - rgrid[1]

        funhr = [h_weight(r, angs, sl, da, taper_a) for r in rgrid]

        return [2π * sum( besselj0(2*pi*ti*rgrid[i]) * rgrid[i] * funhr[i] for i in eachindex(rgrid)) * dr for ti in t]
    end

    """
        h_weight(r, angs, sl, da, taper_a)

    helper function for taper bias
    """
    function h_weight(r, angs, sl, da, taper_a)
        # taper-weighted isotropized set_cov i.e. average of h*set_cov over directions.
        u = [(r*cos(a), r*sin(a)) for a in angs]
        set_cov = [max((sl[1]-abs(ui[1])) * (sl[2]-abs(ui[2])), 0) for ui in u]
        2sum(sqexp_taper(u[i], sl, taper_a) * set_cov[i] for i in eachindex(u)) * da / (2*pi)
    end

    ## taper
    """
        sqexp_taper(x, sl, taper_a::Real)

    isotropic taper
    """
    function sqexp_taper(x, sl, taper_a::Real)
        V = prod(sl)
        exp( -taper_a/4 * ((x[1]/sl[1])^2 + (x[2]/sl[2])^2 )  ) / V
    end

    
    ## internal sum
    """
        iso_sum(x::PointSet, t::Real, sl, taper_a::Real)

    Computes the isotropic transform for a single choice of freq radius t
    """
    function iso_sum(x::PointSet, t::AbstractVector, sl, taper_a::Real)
        pts = x.geoms

        out = zeros(length(t))
        for i in firstindex(pts):lastindex(pts)-1, j in i+1:lastindex(pts)
            Δ = pts[i].coords .- pts[j].coords
            d = sqrt(Δ[1]^2+Δ[2]^2) * 2π
            hij = sqexp_taper(Δ, sl, taper_a)
            for k in eachindex(out)
                out[k] += besselj0(d * t[k]) * hij
            end
        end

        return out
    end
end