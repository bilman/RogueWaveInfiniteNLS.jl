# module rwio
# rwio.jl
# include("rwio_RHP.jl")
# using rwio_RHP

# The changes below are due to the renaming the number of points used on each arc in all the rwio_foo() routines.
# global const NODEF_PTS = 800
global const NODEF_PTS = 400
# global const LARGEX_PTS = 280
global const LARGEX_PTS = 140
# global const LARGET_PTS = 280
global const LARGET_PTS = 140

# Compute RWIO Q(X,T;a,b,\beta) in the first quadrant Q1

function psi_undeformed_Q1(X,T,a,b,β,npts::Integer)
    # Xsc = β*X
    # Tsc = β^2*T
    return(β*rwio_undeformed(β*X,β^2*T,a,b,npts))
end

function psi_largeX_Q1(X,T,a,b,β,npts::Integer)
    # Xsc = β*X
    # Tsc = β^2*T
    return β*rwio_largeX(β*X,vfromXT(β*X,β^2*T),a,b, npts::Integer)
end

function psi_Painleve_Q1(X,T,a,b,β,npts::Integer)
    # Xsc = β*X
    # Tsc = β^2*T
    # v = vfXT(Xsc,Tsc)
    return β*rwio_Painleve(β*X,vfromXT(β*X,β^2*T),a,b,npts)
end

function psi_largeT_Q1(X,T,a,b,β, npts::Integer)
    # Xsc = β*X
    # Tsc = β^2*T
    return β*rwio_largeT(β^2*T,wfromXT(β*X,β^2*T),a,b,npts)
end

# Compute RWIO Q(X,T;a,b,\beta) in the second quadrant Q2

function psi_undeformed_Q2(X,T,a,b,β,npts::Integer)
    # Xsc = β*X
    # Tsc = β^2*T
    return(β*rwio_undeformed(β*abs(X),β^2*T,b,a,npts))
end

function psi_largeX_Q2(X,T,a,b,β,npts::Integer)
    # Xsc = β*X
    # Tsc = β^2*T
    return β*rwio_largeX(β*abs(X),vfromXT(β*abs(X),β^2*T),b,a, npts::Integer)
end

function psi_Painleve_Q2(X,T,a,b,β,npts::Integer)
    # Xsc = β*X
    # Tsc = β^2*T
    # v = vfXT(Xsc,Tsc)
    return β*rwio_Painleve(β*abs(X),vfromXT(β*abs(X),β^2*T),b,a,npts)
end

function psi_largeT_Q2(X,T,a,b,β, npts::Integer)
    # Xsc = β*X
    # Tsc = β^2*T
    return β*rwio_largeT(β^2*T,wfromXT(β*abs(X),β^2*T),b,a,npts)
end

# Compute RWIO Q(X,T;a,b,\beta) in the fourth quadrant Q4

function psi_undeformed_Q4(X,T,a,b,β, npts::Integer)
    # Xsc = β*X
    # Tsc = β^2*T
    return conj(β*rwio_undeformed(β*X,β^2*abs(T),conj(a),conj(b),npts))
end

function psi_largeX_Q4(X,T,a,b,β,npts::Integer)
    # Xsc = β*X
    # Tsc = β^2*T
    return conj(β*rwio_largeX(β*X,vfromXT(β*X,β^2*abs(T)),conj(a),conj(b),npts))
end

function psi_Painleve_Q4(X,T,a,b,β,npts::Integer)
    # Xsc = β*X
    # Tsc = β^2*T
    return conj(β*rwio_Painleve(β*X,vfromXT(β*X,β^2*abs(T)),conj(a),conj(b),npts))
end

function psi_largeT_Q4(X,T,a,b,β,npts::Integer)
    return conj(β*rwio_largeT(β^2*abs(T),wfromXT(β*X,β^2*abs(T)),conj(a),conj(b),npts))
end

# Compute RWIO Q(X,T;a,b,\beta) in the third quadrant Q3

function psi_undeformed_Q3(X,T,a,b,β,numpts::Integer)
    # Xsc = β*X
    # Tsc = β^2*T
    return conj(β*rwio_undeformed(β*abs(X),β^2*abs(T),conj(b),conj(a),numpts))
end

function psi_largeX_Q3(X,T,a,b,β,numpts::Integer)
    # Xsc = β*X
    # Tsc = β^2*T
    return conj(β*rwio_largeX(β*abs(X),vfromXT(β*abs(X),β^2*abs(T)),conj(b),conj(a),numpts))
end

function psi_Painleve_Q3(X,T,a,b,β,numpts::Integer)
    # Xsc = β*X
    # Tsc = β^2*T
    return conj(β*rwio_Painleve(β*abs(X),vfromXT(β*abs(X),β^2*abs(T)),conj(b),conj(a),numpts))
end

function psi_largeT_Q3(X,T,a,b,β,numpts::Integer)
    # Xsc = β*X
    # Tsc = β^2*T
    return conj(β*rwio_largeT(β^2*abs(T),wfromXT(β*abs(X),β^2*abs(T)),conj(b),conj(a),numpts))
end

# Combine the routines in the quadrants to compute Q(X,T,[a,b],β)

function psi_undeformed(X,T,a,b,β,numpts::Integer)
    # Xsc = β*X
    # Tsc = β^2*T
    # return β*nlsnodeformationRescaledT(β*X,β^2*T,a,b)
    # Scaling is done in the routines that compute for the given quadrants.
    # So we just pass (X,T,β) as they are: Given (X,T) and some β>0.
    if T>=0 && X>=0
        return psi_undeformed_Q1(X,T,a,b,β, numpts)
    elseif T>=0 && X<0
        return psi_undeformed_Q2(X,T,a,b,β, numpts)
    elseif T<0 && X<0
        return psi_undeformed_Q3(X,T,a,b,β, numpts)
    elseif T<0 && X>=0
        return psi_undeformed_Q4(X,T,a,b,β, numpts)
    end
end

function psi_largeX(X,T,a,b,β,numpts::Integer)
    # Scaling is done in the routines that compute for the given quadrants.
    # So we just pass (X,T,β) as they are: Given (X,T) and some β>0.
    if T>=0 && X>=0
        return psi_largeX_Q1(X,T,a,b,β,numpts)
    elseif T>=0 && X<0
        return psi_largeX_Q2(X,T,a,b,β,numpts)
    elseif T<0 && X<0
        return psi_largeX_Q3(X,T,a,b,β,numpts)
    elseif T<0 && X>=0
        return psi_largeX_Q4(X,T,a,b,β,numpts)
    end
end


function psi_Painleve(X,T,a,b,β,numpts::Integer)
    # Scaling is done in the routines that compute for the given quadrants.
    # So we just pass (X,T,β) as they are: Given (X,T) and some β>0.
    if T>=0 && X>=0
        return psi_Painleve_Q1(X,T,a,b,β,numpts)
    elseif T>=0 && X<0
        return psi_Painleve_Q2(X,T,a,b,β,numpts)
    elseif T<0 && X<0
        return psi_Painleve_Q3(X,T,a,b,β,numpts)
    elseif T<0 && X>=0
        return psi_Painleve_Q4(X,T,a,b,β,numpts)
    end
end

function psi_largeT(X,T,a,b,β,numpts::Integer)
    # Scaling is done in the routines that compute for the given quadrants.
    # So we just pass (X,T,β) as they are: Given (X,T) and some β>0.
    if T>=0 && X>=0
        return psi_largeT_Q1(X,T,a,b,β,numpts)
    elseif T>=0 && X<0
        return psi_largeT_Q2(X,T,a,b,β,numpts)
    elseif T<0 && X<0
        return psi_largeT_Q3(X,T,a,b,β,numpts)
    elseif T<0 && X>=0
        return psi_largeT_Q4(X,T,a,b,β,numpts)
    end
end

function psi(X,T,a,b,β)
    if T>125000
        println("T is too large. Accuracy will be lost.")
        return()
    end
    # This implements contour rescaling (without changing the solution) in the
    # no-deformation setting for T<=8. We divide the circle jump
    # by sqrt(T) so that the term T*λ^2 in the phase behaves as if T=1.
    Xsc = β*X
    Tsc = β^2*T
    v = vfromXT(abs(Xsc), abs(Tsc))
    w = wfromXT(abs(Xsc), abs(Tsc))
    # Never pass in the rescaled X,T values. The underlying routines work with those.
    if Xsc^2 + Tsc^2 <= 4.
        # (X,T) is in a small disk, so no deformation needed.
        # println("No deformation")
        return psi_undeformed(X, T, a, b, β, NODEF_PTS)
    elseif abs(Tsc)<=8.
        # Below the level T=8. No deformation is used in the large-T domain.
        # But the large-X deformation is employed since this region is a horizontal strip.
        if (v>VCRIT) && (abs(v-VCRIT)>VDISTANCEPAINLEVE)
            # println("Below T=8 but no deformation.")
            psi_undeformed(X, T, a, b, β, NODEF_PTS)
        elseif (v<VCRIT) && (abs(v-VCRIT)>VDISTANCEPAINLEVE)
            # println("Below T=8 and large-X deformation.")
            psi_largeX(X, T, a, b, β, LARGEX_PTS)
        else
            # println("Below T=8 and Painlevé deformation.")
            psi_Painleve(X, T, a, b, β, LARGEX_PTS)
        end
    else
        # Above the level T=8. All deformations are employed
        # But the large-X deformation is employed.
        # if (v>VCRIT) && (abs(v-VCRIT)>VDISTANCEPAINLEVE)
        if (v>VCRIT) && (w<0.98*WCRIT)
            # println("Above T=8 and large-T deformation.")
            psi_largeT(X, T, a, b, β, LARGET_PTS)
        elseif (v<VCRIT) && (abs(v-VCRIT)>VDISTANCEPAINLEVE)
            # println("Above T=8 and large-X deformation.")
            psi_largeX(X, T, a, b, β, LARGEX_PTS)
        else
            # println("Above T=8 and Painlevé deformation.")
            psi_Painleve(X, T, a, b, β, LARGEX_PTS)
        end
    end
end

# This is used for T=0 to compute the elliptic spread.
function Psi_slice_T0(X,a,b,β)
    if X<=8.
        return psi_undeformed(X,0,a,b,β,NODEF_PTS)
    else
        return psi_largeX(X,0,a,b,β,LARGEX_PTS)
    end
end



# end # module
