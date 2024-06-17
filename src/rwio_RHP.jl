# module rwio_RHP
# rwio_RHP.jl

# export wfromXT, vfromXT, XfromTw, vfromTw, Tmax_givenX, Xmax_givenT, rwio_nodeformation_RHP, rwio_nodeformation, rwio_largeX_RHP, rwio_largeX, rwio_Painleve_RHP, rwio_Painleve, rwio_largeT_RHP, rwio_largeT, WCRIT, VCRIT, VDISTANCEPAINLEVE

using OperatorApproximation
using LinearAlgebra
using Polynomials

global const WCRIT = 54^(1/3)
global const VCRIT = 1/sqrt(54)
global const WDISTANCEPAINLEVE = 0.0185
global const VDISTANCEPAINLEVE = 0.00025
global const P2SWITCHCONST = 1
global const VKILL1 = VCRIT*0.9
global const VKILL2 = VCRIT-VDISTANCEPAINLEVE


# Plotting routines

function rhwellposed_show(rhp)
    l = length(rhwellposed(rhp))
    empty_bool_array = zeros(Bool, l)
    for i in 1:l
        empty_bool_array[i] = rhwellposed(rhp)[i][2]≈[1 0; 0 1]
    end
    return empty_bool_array
end

# General routines

function wfromXT(X,T)
    return X/T^(2/3)
end

function vfromXT(X,T)
    return T/X^(3/2)
end

function XfromTw(T,w)
    return w*T^(2/3)
end

function TfromXv(X,v)
    return X^(3/2)*v
end

function vfromTw(T,w) 
    return T/(XfromTw(T,w))^(3/2)
end

function Tmax_givenX(X)
    return X^(3/2)*VCRIT
end

function Xmax_givenT(T)
    return T^(2/3)*WCRIT
end

# function painleve_switch(X,T)
#     v = vfromXT(X,T)
#     if X>10
#     else
#     end
#     return(abs(vfromXT(X,T)))
# end

function arcpoints_largeX(X,v)
    # This function never gets called (and shouldn't) if VCRIT-v>=VDISTANCEPAINLEVE
    # It returns the vertices of the polygon in the UHP that models the jump contour outside the circle.
    # θ = z -> z + v*z^2 + 2/z
    if 0< v && v<= VCRIT/4
        z1 =  (1 / (6 * v)) * (-1 + 2 * cos((-1/3) * asin(v * sqrt(216 * (1 - 54 * v^2))) - (1/3) * π))
        z2 =  (1 / (6 * v)) * (-1 + 2 * cos((-1/3) * asin(v * sqrt(216 * (1 - 54 * v^2))) + (1/3) * π))
    elseif VCRIT/4<v && v < VCRIT
        z1 = (1 / (6 * v)) * (-1 + 2 * cos((1/3) * acos(108 * v^2 - 1) - (2/3) * π))
        z2 = (1 / (6 * v)) * (-1 + 2 * cos((1/3) * acos(108 * v^2 - 1) ))
    elseif v==0
        z1 = -sqrt(2)
        z2 = sqrt(2)
    else
        println("Bad v.")
        return
    end
    outrad = 2
    inrad = 1
    if 0<=v<VKILL1
        CONFIG = 1  # This is a flag that gets returned.
        # zinf=rootlist[3] is not needed.
        arcpointsOut = [z1, z1 + outrad*exp(3im*pi/4), 0.5*(z1+z2) + outrad*1im, z2 + outrad*exp(1im*pi/4), z2]
    elseif VKILL1<=v<min(VKILL2, X^(-1/3))
        CONFIG = 2
        polyθ = Polynomial([-2,0,1,2*v])
        rootlist = Polynomials.roots(polyθ) |> reverse
        zinf=rootlist[3]
        z1leftrad = min(0.6*abs(zinf-z1)*sqrt(2),outrad)
        arcpointsOut = [z1, z1 + z1leftrad*exp(3im*pi/4), real(z1 + z1leftrad*exp(3im*pi/4))+outrad*1im, 0.5*(z1+z2) + outrad*1im, z2 + outrad*exp(1im*pi/4), z2]
    else
        println("Bad v value for arc determination!")
        CONFIG = 0
        return 
    end
    arcpointsIn = [z1, z1 + inrad*exp(1im*pi/4), 0.5*(z1+z2) + inrad*1im, z2 + inrad*exp(3im*pi/4), z2]
    return(arcpointsOut, arcpointsIn, CONFIG)
end

function angle_z0z1(w)
    # This function computes the argument (phase, direction) of the local path
    # that coincides with the branch cut of the g-function in the paper
    # emanating from z0 towards z1. This phase factor then gets used to rotate the polygon at z0
    # in rwio_largeT.
    aw = w/6.
    b = -(2.0^(-1/3))
    cw = 2*w/3.0
    d = (2.0)^(2/3)
    bp1w = 0.5*(-cw + 1im*sqrt(4*d - cw^2))
    bp2w = 0.5*(-cw - 1im*sqrt(4*d - cw^2))
    hprime = z -> (2/z^2)*(z^2 + aw*z + b)*(z - bp2w)*sqrt((z - bp1w)/(z - bp2w))
    ds = 0.0001
    arclen = 0.
    direction = 0. + 1im
    pts = complex(zeros(0))
    append!(pts, 0.5*(-aw - sqrt(aw^2 - 4*b)))
    append!(pts, pts[end] + direction*ds)
    while (abs(pts[end]-bp1w)>100*ds && arclen<4)
        direction = -abs(hprime(pts[end]))/hprime(pts[end])
        append!(pts, pts[end] + direction*ds)
        arclen = arclen + ds
    end
    return(angle(bp1w - pts[end]))
end

# The RH solver operates on matrices of functions rather than matrix-valued functions.
# Below is a built-in routine that takes a matrix-valued function and returns a 2x2 matrix of scalar-valued functions.
mkfun = y -> mvf2mof(y,2,2)


# Compute RWIO without any deformation

function rwio_nodeformation_RHP(X,T,a,b,npts::Integer)
    # tau = abs(b/conj(a))
    N = sqrt(abs(a)^2+abs(b)^2)
    θ = z->X*z+T*z^2+2/z

    G11 = z->a/N
    G12 = z->(conj(b)/N)*exp(-2im*θ(z))
    G21 = z->(-b/N)*exp(2im*θ(z))
    G22 = z->conj(a)/N

    G = [G11 G12; G21 G22]
    # Unit circle modeled by a square
    intervals = [exp(1im*pi/4) exp(-1im*pi/4); exp(-1im*pi/4) exp(-3im*pi/4); exp(-3im*pi/4) exp(3im*pi/4); exp(3im*pi/4) exp(1im*pi/4)]
    jumps = [G, G, G, G]
    rhp = RHP(intervals, jumps)
    rhp_solver = RHSolver(rhp)
    numarcs = length(jumps)
    u1 = rhp_solver([1 0],numarcs*npts)
    u2 = rhp_solver([0 1],numarcs*npts)
    U = rhp_solver(([1 0], [0 1]), numarcs*npts)
    U = hcat(U...) |> transpose |> Matrix
    Usums = sum.(U)
    nls = -(1/π)*Usums[1,2]
    return nls, rhp, rhp_solver, U, u1, u2
end


function rwio_nodeformation(X,T,a,b,npts::Integer)
    # tau = abs(b/conj(a))
    N = sqrt(abs(a)^2+abs(b)^2)
    θ = z->X*z+T*z^2+2/z

    G11 = z->a/N
    G12 = z->(conj(b)/N)*exp(-2im*θ(z))
    G21 = z->(-b/N)*exp(2im*θ(z))
    G22 = z->conj(a)/N

    G = [G11 G12; G21 G22]
    # Unit circle modeled by a square
    intervals = [exp(1im*pi/4) exp(-1im*pi/4); exp(-1im*pi/4) exp(-3im*pi/4); exp(-3im*pi/4) exp(3im*pi/4); exp(3im*pi/4) exp(1im*pi/4)]
    jumps = [G, G, G, G]
    rhp = RHP(intervals, jumps)
    rhp_solver = RHSolver(rhp)
    numarcs = length(jumps)
    # u1 = rhp_solver([1 0],numarcs*npts)
    # u2 = rhp_solver([0 1],numarcs*npts)
    U = rhp_solver(([1 0], [0 1]), numarcs*npts)
    GC.gc()
    U = hcat(U...) |> transpose |> Matrix
    Usums = sum.(U)
    nls = -(1/π)*Usums[1,2];
    return nls
end

function rwio_nodeformation_rescaled_RHP(X,T,a,b,npts::Integer)
    # tau = abs(b/conj(a))
    N = sqrt(abs(a)^2+abs(b)^2)
    θ = z->X*z+T*z^2+2/z
    if abs(T)<=1.
        rad = 1
    else
        rad = 1/sqrt(abs(T))
    end

    G11 = z->a/N
    G12 = z->(conj(b)/N)*exp(-2im*θ(z))
    G21 = z->(-b/N)*exp(2im*θ(z))
    G22 = z->conj(a)/N

    G = [G11 G12; G21 G22]
    # Unit circle modeled by a square
    intervals = [exp(1im*pi/4)*rad exp(-1im*pi/4)*rad; exp(-1im*pi/4)*rad exp(-3im*pi/4)*rad; exp(-3im*pi/4)*rad exp(3im*pi/4)*rad; exp(3im*pi/4)*rad exp(1im*pi/4)*rad]
    jumps = [G, G, G, G]
    rhp = RHP(intervals, jumps)
    rhp_solver = RHSolver(rhp)
    numarcs = length(jumps)
    u1 = rhp_solver([1 0],numarcs*npts)
    u2 = rhp_solver([0 1],numarcs*npts)
    U = rhp_solver(([1 0], [0 1]), numarcs*npts)
    U = hcat(U...) |> transpose |> Matrix
    Usums = sum.(U)
    nls = -(1/π)*Usums[1,2]
    return nls, rhp, rhp_solver, U, u1, u2
    # return Usums, Usums
end

function rwio_nodeformation_rescaled(X,T,a,b,npts::Integer)
    # tau = abs(b/conj(a))
    N = sqrt(abs(a)^2+abs(b)^2)
    θ = z->X*z+T*z^2+2/z
    if abs(T)<=1.
        rad = 1
    else
        rad = 1/sqrt(abs(T))
    end

    G11 = z->a/N
    G12 = z->(conj(b)/N)*exp(-2im*θ(z))
    G21 = z->(-b/N)*exp(2im*θ(z))
    G22 = z->conj(a)/N

    G = [G11 G12; G21 G22]
    # Unit circle modeled by a square
    intervals = [exp(1im*pi/4)*rad exp(-1im*pi/4)*rad; exp(-1im*pi/4)*rad exp(-3im*pi/4)*rad; exp(-3im*pi/4)*rad exp(3im*pi/4)*rad; exp(3im*pi/4)*rad exp(1im*pi/4)*rad]
    jumps = [G, G, G, G]
    rhp = RHP(intervals, jumps)
    rhp_solver = RHSolver(rhp)
    numarcs = length(jumps)
    # u1 = rhp_solver([1 0],numarcs*npts)
    # u2 = rhp_solver([0 1],numarcs*npts)
    U = rhp_solver(([1 0], [0 1]), numarcs*npts)
    GC.gc()
    U = hcat(U...) |> transpose |> Matrix
    Usums = sum.(U)
    nls = -(1/π)*Usums[1,2]
    return nls
    # return Usums, Usums
end

# Compute RWIO with large-X deformations

function rwio_largeX_RHP(X,v,a,b, npts::Integer)
    XX = sqrt(X)
    if v >= VCRIT
        println("v value is too large.")
        return
    end
    afrak = abs(a)/sqrt(abs(a)^2 + abs(b)^2)
    bfrak = abs(b)/sqrt(abs(a)^2 + abs(b)^2)
    τ = bfrak/afrak
    θ = z->z+v*z^2+2/z
    if 0< v && v<= VCRIT/4
        z1 =  (1 / (6 * v)) * (-1 + 2 * cos((-1/3) * asin(v * sqrt(216 * (1 - 54 * v^2))) - (1/3) * π))
        z2 =  (1 / (6 * v)) * (-1 + 2 * cos((-1/3) * asin(v * sqrt(216 * (1 - 54 * v^2))) + (1/3) * π))
    elseif VCRIT/4<v && v < VCRIT
        z1 = (1 / (6 * v)) * (-1 + 2 * cos((1/3) * acos(108 * v^2 - 1) - (2/3) * π))
        z2 = (1 / (6 * v)) * (-1 + 2 * cos((1/3) * acos(108 * v^2 - 1) ))
    elseif v==0
        z1 = -sqrt(2)
        z2 = sqrt(2)
    else
        println("Bad value of v.")
        return
    end
    Lup = [z->1 z->0; z->-τ*exp(2im*XX*θ(z)) z->1]
    Rup = [z->1 z->(bfrak*afrak)*exp(-2im*XX*θ(z)); z->0 z->1]
    Djump = [z->afrak^(-2) z->0; z->0 z->afrak^2]
    Rdown = [z->1 z->0; z->-(bfrak*afrak)*exp(2im*XX*θ(z)) z->1]
    Ldown = [z->1 z->τ*exp(-2im*XX*θ(z)); z->0 z->1]

    outptsUHP, inptsUHP, CONFIG = arcpoints_largeX(X,v)
    outptsLHP = outptsUHP|>reverse|>conj
    inptsLHP = inptsUHP|>reverse|>conj
    if CONFIG == 1
        arcsUHPout = [outptsUHP[1] outptsUHP[2]; outptsUHP[2] outptsUHP[3]; outptsUHP[3] outptsUHP[4]; outptsUHP[4] outptsUHP[5]]
        arcsLHPout = [outptsLHP[1] outptsLHP[2]; outptsLHP[2] outptsLHP[3]; outptsLHP[3] outptsLHP[4]; outptsLHP[4] outptsLHP[5]]
        J = [Lup,Lup,Lup,Lup,Ldown,Ldown,Ldown,Ldown,Rup,Rup,Rup,Rup,Rdown,Rdown,Rdown,Rdown,Djump]
    elseif CONFIG == 2
        arcsUHPout = [outptsUHP[1] outptsUHP[2]; outptsUHP[2] outptsUHP[3]; outptsUHP[3] outptsUHP[4]; outptsUHP[4] outptsUHP[5]; outptsUHP[5] outptsUHP[6]]
        arcsLHPout = [outptsLHP[1] outptsLHP[2]; outptsLHP[2] outptsLHP[3]; outptsLHP[3] outptsLHP[4]; outptsLHP[4] outptsLHP[5]; outptsLHP[5] outptsLHP[6]]
        J = [Lup,Lup,Lup,Lup,Lup,Ldown,Ldown,Ldown,Ldown,Ldown,Rup,Rup,Rup,Rup,Rdown,Rdown,Rdown,Rdown,Djump]
    else
        println("Problem with CONFIG in large-X arc setup.")
        return
    end
    arcsUHPin = [inptsUHP[1] inptsUHP[2]; inptsUHP[2] inptsUHP[3]; inptsUHP[3] inptsUHP[4]; inptsUHP[4] inptsUHP[5]]
    arcsLHPin = [inptsLHP[1] inptsLHP[2]; inptsLHP[2] inptsLHP[3]; inptsLHP[3] inptsLHP[4]; inptsLHP[4] inptsLHP[5]]
    arcReal = [inptsUHP[1] inptsUHP[5]]

    intervals = vcat(arcsUHPout, arcsLHPout, arcsUHPin, arcsLHPin, arcReal)
    rhp = RHP(intervals, J)
    rhp_solver = RHSolver(rhp)
    numarcs = length(J)
    u1 = rhp_solver([1 0],numarcs*npts)
    u2 = rhp_solver([0 1],numarcs*npts)
    U = rhp_solver(([1 0], [0 1]), numarcs*npts)
    U = hcat(U...) |> transpose |> Matrix
    Usums = sum.(U)
    nls = -(1/π)*Usums[1,2]
    nls = nls*exp(-1im*angle(a*b))/XX
    return nls, rhp, rhp_solver, U, u1, u2
end


function rwio_largeX(X,v,a,b, npts::Integer)
    XX = sqrt(X)
    if v >= VCRIT
        println("v value is too large.")
        return
    end
    afrak = abs(a)/sqrt(abs(a)^2 + abs(b)^2)
    bfrak = abs(b)/sqrt(abs(a)^2 + abs(b)^2)
    τ = bfrak/afrak
    θ = z->z+v*z^2+2/z
    if 0< v && v<= VCRIT/4
        z1 =  (1 / (6 * v)) * (-1 + 2 * cos((-1/3) * asin(v * sqrt(216 * (1 - 54 * v^2))) - (1/3) * π))
        z2 =  (1 / (6 * v)) * (-1 + 2 * cos((-1/3) * asin(v * sqrt(216 * (1 - 54 * v^2))) + (1/3) * π))
    elseif VCRIT/4<v && v < VCRIT
        z1 = (1 / (6 * v)) * (-1 + 2 * cos((1/3) * acos(108 * v^2 - 1) - (2/3) * π))
        z2 = (1 / (6 * v)) * (-1 + 2 * cos((1/3) * acos(108 * v^2 - 1) ))
    elseif v==0
        z1 = -sqrt(2)
        z2 = sqrt(2)
    else
        println("Bad value of v.")
        return
    end
    Lup = [z->1 z->0; z->-τ*exp(2im*XX*θ(z)) z->1]
    Rup = [z->1 z->(bfrak*afrak)*exp(-2im*XX*θ(z)); z->0 z->1]
    Djump = [z->afrak^(-2) z->0; z->0 z->afrak^2]
    Rdown = [z->1 z->0; z->-(bfrak*afrak)*exp(2im*XX*θ(z)) z->1]
    Ldown = [z->1 z->τ*exp(-2im*XX*θ(z)); z->0 z->1]

    outptsUHP, inptsUHP, CONFIG = arcpoints_largeX(X,v)
    outptsLHP = outptsUHP|>reverse|>conj
    inptsLHP = inptsUHP|>reverse|>conj
    if CONFIG == 1
        arcsUHPout = [outptsUHP[1] outptsUHP[2]; outptsUHP[2] outptsUHP[3]; outptsUHP[3] outptsUHP[4]; outptsUHP[4] outptsUHP[5]]
        arcsLHPout = [outptsLHP[1] outptsLHP[2]; outptsLHP[2] outptsLHP[3]; outptsLHP[3] outptsLHP[4]; outptsLHP[4] outptsLHP[5]]
        J = [Lup,Lup,Lup,Lup,Ldown,Ldown,Ldown,Ldown,Rup,Rup,Rup,Rup,Rdown,Rdown,Rdown,Rdown,Djump]
    elseif CONFIG == 2
        arcsUHPout = [outptsUHP[1] outptsUHP[2]; outptsUHP[2] outptsUHP[3]; outptsUHP[3] outptsUHP[4]; outptsUHP[4] outptsUHP[5]; outptsUHP[5] outptsUHP[6]]
        arcsLHPout = [outptsLHP[1] outptsLHP[2]; outptsLHP[2] outptsLHP[3]; outptsLHP[3] outptsLHP[4]; outptsLHP[4] outptsLHP[5]; outptsLHP[5] outptsLHP[6]]
        J = [Lup,Lup,Lup,Lup,Lup,Ldown,Ldown,Ldown,Ldown,Ldown,Rup,Rup,Rup,Rup,Rdown,Rdown,Rdown,Rdown,Djump]
    else
        println("Problem with CONFIG in large-X arc setup.")
        return
    end
    arcsUHPin = [inptsUHP[1] inptsUHP[2]; inptsUHP[2] inptsUHP[3]; inptsUHP[3] inptsUHP[4]; inptsUHP[4] inptsUHP[5]]
    arcsLHPin = [inptsLHP[1] inptsLHP[2]; inptsLHP[2] inptsLHP[3]; inptsLHP[3] inptsLHP[4]; inptsLHP[4] inptsLHP[5]]
    arcReal = [inptsUHP[1] inptsUHP[5]]

    intervals = vcat(arcsUHPout, arcsLHPout, arcsUHPin, arcsLHPin, arcReal)
    rhp = RHP(intervals, J)
    rhp_solver = RHSolver(rhp)
    numarcs = length(J)
    # u1 = rhp_solver([1 0],numarcs*npts)
    # u2 = rhp_solver([0 1],numarcs*npts)
    U = rhp_solver(([1 0], [0 1]), numarcs*npts)
    U = hcat(U...) |> transpose |> Matrix
    Usums = sum.(U)
    nls = -(1/π)*Usums[1,2];
    nls = nls*exp(-1im*angle(a*b))/XX
    return nls
end

function rwio_largeX_NOFRAKS(X,v,a,b, npts::Integer)
    XX = sqrt(X)
    if v >= VCRIT
        println("v value is too large.")
        return
    end
    τ = abs(b/a)
    ν = angle(b/conj(a))
    p = (1/(2*pi))*log(1+τ^2)
    θ = z->z+v*z^2+2/z
    if 0< v && v<= VCRIT/4
        z1 =  (1 / (6 * v)) * (-1 + 2 * cos((-1/3) * asin(v * sqrt(216 * (1 - 54 * v^2))) - (1/3) * π))
        z2 =  (1 / (6 * v)) * (-1 + 2 * cos((-1/3) * asin(v * sqrt(216 * (1 - 54 * v^2))) + (1/3) * π))
    elseif VCRIT/4<v && v < VCRIT
        z1 = (1 / (6 * v)) * (-1 + 2 * cos((1/3) * acos(108 * v^2 - 1) - (2/3) * π))
        z2 = (1 / (6 * v)) * (-1 + 2 * cos((1/3) * acos(108 * v^2 - 1) ))
    elseif v==0
        z1 = -sqrt(2)
        z2 = sqrt(2)
    else
        println("Bad value for v.")
        return
    end
    Lup = [z->1 z->0; z->-τ*exp(2im*XX*θ(z)) z->1]
    Rup = [z->1 z->τ*exp(-2*pi*p)*exp(-2im*XX*θ(z)); z->0 z->1]
    Djump = [z->(1+τ^2) z->0; z->0 z->1/(1+τ^2)]
    Rdown = [z->1 z->0; z->-τ*exp(-2*pi*p)*exp(2im*XX*θ(z)) z->1]
    Ldown = [z->1 z->τ*exp(-2im*XX*θ(z)); z->0 z->1]

    outptsUHP, inptsUHP, CONFIG = arcpoints_largeX(X,v)
    outptsLHP = outptsUHP|>reverse|>conj
    inptsLHP = inptsUHP|>reverse|>conj
    if CONFIG == 1
        arcsUHPout = [outptsUHP[1] outptsUHP[2]; outptsUHP[2] outptsUHP[3]; outptsUHP[3] outptsUHP[4]; outptsUHP[4] outptsUHP[5]]
        arcsLHPout = [outptsLHP[1] outptsLHP[2]; outptsLHP[2] outptsLHP[3]; outptsLHP[3] outptsLHP[4]; outptsLHP[4] outptsLHP[5]]
        J = [Lup,Lup,Lup,Lup,Ldown,Ldown,Ldown,Ldown,Rup,Rup,Rup,Rup,Rdown,Rdown,Rdown,Rdown,Djump]
    elseif CONFIG == 2
        arcsUHPout = [outptsUHP[1] outptsUHP[2]; outptsUHP[2] outptsUHP[3]; outptsUHP[3] outptsUHP[4]; outptsUHP[4] outptsUHP[5]; outptsUHP[5] outptsUHP[6]]
        arcsLHPout = [outptsLHP[1] outptsLHP[2]; outptsLHP[2] outptsLHP[3]; outptsLHP[3] outptsLHP[4]; outptsLHP[4] outptsLHP[5]; outptsLHP[5] outptsLHP[6]]
        J = [Lup,Lup,Lup,Lup,Lup,Ldown,Ldown,Ldown,Ldown,Ldown,Rup,Rup,Rup,Rup,Rdown,Rdown,Rdown,Rdown,Djump]
    else
        println("Problem with CONFIG in large-X arc setup.")
        return
    end
    arcsUHPin = [inptsUHP[1] inptsUHP[2]; inptsUHP[2] inptsUHP[3]; inptsUHP[3] inptsUHP[4]; inptsUHP[4] inptsUHP[5]]
    arcsLHPin = [inptsLHP[1] inptsLHP[2]; inptsLHP[2] inptsLHP[3]; inptsLHP[3] inptsLHP[4]; inptsLHP[4] inptsLHP[5]]
    arcReal = [inptsUHP[1] inptsUHP[5]]

    intervals = vcat(arcsUHPout, arcsLHPout, arcsUHPin, arcsLHPin, arcReal)
    rhp = RHP(intervals, J)
    rhp_solver = RHSolver(rhp)
    numarcs = length(J)
    # u1 = rhp_solver([1 0],numarcs*npts)
    # u2 = rhp_solver([0 1],numarcs*npts)
    U = rhp_solver(([1 0], [0 1]), numarcs*npts)
    GC.gc()
    U = hcat(U...) |> transpose |> Matrix
    Usums = sum.(U)
    nls = -(1/π)*Usums[1,2];
    nls = nls*exp(-1im*ν)/XX
    return nls
end


function rwio_Painleve_RHP(X,v,a,b,npts::Integer)
    XX = sqrt(X)
    afrak = abs(a)/sqrt(abs(a)^2 + abs(b)^2)
    bfrak = abs(b)/sqrt(abs(a)^2 + abs(b)^2)
    τ = bfrak/afrak
    θ = z->z+v*z^2+2/z
    # T = v*X^(3/2)
    z1 = -sqrt(6)
    z2 = (1 / (6 * v)) * (-1 + 2 * cos((1/3) * acos((108 * v^2 - 1)|>ComplexF64) )) |> real
    outrad = 2
    inrad = 1
    # Jumps on lenses
    Lup = [z->1 z->0; z->-τ*exp(2im*XX*θ(z)) z->1]
    Rup = [z->1 z->(bfrak*afrak)*exp(-2im*XX*θ(z)); z->0 z->1]
    Djump = [z->afrak^(-2) z->0; z->0 z->afrak^2]
    Rdown = [z->1 z->0; z->-(bfrak*afrak)*exp(2im*XX*θ(z)) z->1]
    Ldown = [z->1 z->τ*exp(-2im*XX*θ(z)); z->0 z->1]
    Djump = [z->(1+τ^2) z->0; z->0 z->1/(1+τ^2)]

    # Contour endpoints
    outptsUHP = [z1, z1+outrad*exp(1im*pi/2), z2+outrad*sqrt(2)*exp(1im*pi/4),z2]
    arcsUHPout = [outptsUHP[1] outptsUHP[2]; outptsUHP[2] outptsUHP[3]; outptsUHP[3] outptsUHP[4]]
    outptsLHP = conj.(outptsUHP)|>reverse
    arcsLHPout = [outptsLHP[1] outptsLHP[2]; outptsLHP[2] outptsLHP[3]; outptsLHP[3] outptsLHP[4]]
    inptsUHP = [z1, z1+inrad*exp(1im*pi/6), z2+inrad*exp(3im*pi/4),z2]
    arcsUHPin = [inptsUHP[1] inptsUHP[2]; inptsUHP[2] inptsUHP[3]; inptsUHP[3] inptsUHP[4]]
    inptsLHP = conj.(inptsUHP)|>reverse
    arcsLHPin = [inptsLHP[1] inptsLHP[2]; inptsLHP[2] inptsLHP[3]; inptsLHP[3] inptsLHP[4]]
    arcReal = [inptsUHP[1] inptsUHP[4]]
    # Γ = arcsUHPout ∪ arcsLHPout ∪ arcsUHPin ∪ arcsLHPin ∪ arcReal
    intervals = vcat(arcsUHPout, arcsLHPout, arcsUHPin, arcsLHPin, arcReal)
    J = [Lup,Lup,Lup,Ldown,Ldown,Ldown,Rup,Rup,Rup,Rdown,Rdown,Rdown,Djump]
    rhp = RHP(intervals, J)
    rhp_solver = RHSolver(rhp)
    numarcs = length(J)
    u1 = rhp_solver([1 0],numarcs*npts)
    u2 = rhp_solver([0 1],numarcs*npts)
    U = rhp_solver(([1 0], [0 1]), numarcs*npts)
    U = hcat(U...) |> transpose |> Matrix
    Usums = sum.(U)
    nls = -(1/π)*Usums[1,2]
    nls = nls*exp(-1im*angle(a*b))/XX
    return nls, rhp, rhp_solver, U, u1, u2
end

function rwio_Painleve(X,v,a,b,npts::Integer)
    XX = sqrt(X)
    afrak = abs(a)/sqrt(abs(a)^2 + abs(b)^2)
    bfrak = abs(b)/sqrt(abs(a)^2 + abs(b)^2)
    τ = bfrak/afrak
    θ = z->z+v*z^2+2/z
    # T = v*X^(3/2)
    z1 = -sqrt(6)
    z2 = (1 / (6 * v)) * (-1 + 2 * cos((1/3) * acos((108 * v^2 - 1)|>ComplexF64) )) |> real
    outrad = 2
    inrad = 1
    # Jumps on lenses
    Lup = [z->1 z->0; z->-τ*exp(2im*XX*θ(z)) z->1]
    Rup = [z->1 z->(bfrak*afrak)*exp(-2im*XX*θ(z)); z->0 z->1]
    Djump = [z->afrak^(-2) z->0; z->0 z->afrak^2]
    Rdown = [z->1 z->0; z->-(bfrak*afrak)*exp(2im*XX*θ(z)) z->1]
    Ldown = [z->1 z->τ*exp(-2im*XX*θ(z)); z->0 z->1]
    Djump = [z->(1+τ^2) z->0; z->0 z->1/(1+τ^2)]

    # Contour endpoints
    outptsUHP = [z1, z1+outrad*exp(1im*pi/2), z2+outrad*sqrt(2)*exp(1im*pi/4),z2]
    arcsUHPout = [outptsUHP[1] outptsUHP[2]; outptsUHP[2] outptsUHP[3]; outptsUHP[3] outptsUHP[4]]
    outptsLHP = conj.(outptsUHP)|>reverse
    arcsLHPout = [outptsLHP[1] outptsLHP[2]; outptsLHP[2] outptsLHP[3]; outptsLHP[3] outptsLHP[4]]
    inptsUHP = [z1, z1+inrad*exp(1im*pi/6), z2+inrad*exp(3im*pi/4),z2]
    arcsUHPin = [inptsUHP[1] inptsUHP[2]; inptsUHP[2] inptsUHP[3]; inptsUHP[3] inptsUHP[4]]
    inptsLHP = conj.(inptsUHP)|>reverse
    arcsLHPin = [inptsLHP[1] inptsLHP[2]; inptsLHP[2] inptsLHP[3]; inptsLHP[3] inptsLHP[4]]
    arcReal = [inptsUHP[1] inptsUHP[4]]
    # Γ = arcsUHPout ∪ arcsLHPout ∪ arcsUHPin ∪ arcsLHPin ∪ arcReal
    intervals = vcat(arcsUHPout, arcsLHPout, arcsUHPin, arcsLHPin, arcReal)
    J = [Lup,Lup,Lup,Ldown,Ldown,Ldown,Rup,Rup,Rup,Rdown,Rdown,Rdown,Djump]
    rhp = RHP(intervals, J)
    rhp_solver = RHSolver(rhp)
    numarcs = length(J)
    # u1 = rhp_solver([1 0],numarcs*npts)
    # u2 = rhp_solver([0 1],numarcs*npts)
    U = rhp_solver(([1 0], [0 1]), numarcs*npts)
    GC.gc()
    U = hcat(U...) |> transpose |> Matrix
    Usums = sum.(U)
    nls = -(1/π)*Usums[1,2]
    nls = nls*exp(-1im*angle(a*b))/XX
    return nls
end

# Compute RWIO with large-T deformations

function rwio_largeT_RHP(T,w,a,b, npts::Integer)
    TT = T^(1/3)
    θ = z -> w*z+z^2+2/z
    if w >= WCRIT
        println("w value is too large.")
        return
    end
    af = abs(a)/sqrt(abs(a)^2 + abs(b)^2)
    bf = abs(b)/sqrt(abs(a)^2 + abs(b)^2)
    # τ = bf/af
    Z1 = (1/2)*(-w/6 - sqrt(w^2/36+2^(5/3)))
    Z2 = (1/2)*(-w/6 + sqrt(w^2/36+2^(5/3)))
    Z0 = (1/3)*(-w + 1im*sqrt(WCRIT^2-w^2))
    Z0CC = (1/3)*(-w - 1im*sqrt(WCRIT^2-w^2))
    # argZ0Z1 = angle(Z0-Z1)
    # midZ1Z2 = 0.5*(Z1+Z2)
    # κ is the constant sum of the boundary values of h(z) on Σ.
    κ = -108^(1/3) - (w^2)/3
    # Rotating angle for the polygon at z0.
    argZ0 = angle_z0z1(w)
     # Radii of small disks
    radZ0 = min(1/(1+abs(T)^(2/9)), 0.3*abs(Z1-Z0), 0.25)
    # radZ1 = min(1/(1+abs(T)^(1/6)), 0.3*abs(Z1-Z0), 0.25)

    # Small circle points
    ptsZ0circ = [1,exp(2im*pi/3),exp(1im*pi),exp(-2im*pi/3)]*radZ0*exp(1im*argZ0) .+ Z0
    Z0circ = [ptsZ0circ[1] ptsZ0circ[2]; ptsZ0circ[2] ptsZ0circ[3]; ptsZ0circ[3] ptsZ0circ[4]; ptsZ0circ[4] ptsZ0circ[1]]
    ptsZ0CCcirc = [1,exp(2im*pi/3),exp(1im*pi),exp(-2im*pi/3)]*radZ0*exp(-1im*argZ0) .+ Z0CC
    Z0CCcirc = [ptsZ0CCcirc[1] ptsZ0CCcirc[2]; ptsZ0CCcirc[2] ptsZ0CCcirc[3]; ptsZ0CCcirc[3] ptsZ0CCcirc[4]; ptsZ0CCcirc[4] ptsZ0CCcirc[1]]
    # ptsZ1circ = [1, exp(1im*pi/4), exp(1im*pi/2), exp(3im*pi/4), exp(-3im*pi/4), exp(-1im*pi/2), exp(-1im*pi/4), 1]*radZ1 .+ Z1
    # ptsZ2circ = [exp(1im*pi), exp(-3im*pi/4), exp(-1im*pi/4), exp(1im*pi/4), exp(3im*pi/4)]*radZ1 .+ Z2

    # Internal parameters required for the contours Σ and Γ
    # Exit angle from z0 to the exterior domain
    ϕ0left = angle(ptsZ0circ[2]-Z0)
    # Exit angle from z0 to the interior domain
    ϕ0right = angle(ptsZ0circ[4]-Z0)
    # t1 is the parameter so that Z1+t1*exp(1im*pi/4) intersects the line coming out of Z0 to the interior domain with the Airy angle.
    t1 = (-1/sin(ϕ0right-pi/4))*(imag(Z0)*cos(ϕ0right) + (Z1-real(Z0))*sin(ϕ0right))
    # The midwert below may start with ptsZ1circ instead of Z1 in case there is a little circle at Z1
    # midvert = imag(0.5*(ptsZ1circ[2]+ptsZ0circ[3]))
    # midvert = imag(0.5*(Z1+ptsZ0circ[3]))
    # ΣjunctionInsideUHP = midZ1Z2+1im*midvert
    # ΣjunctionInsideLHP = midZ1Z2-1im*midvert
    ΣjunctionInsideUHP = Z1 + t1*exp(1im*pi/4)
    ΣjunctionInsideLHP = Z1 + t1*exp(-1im*pi/4)

    # Σ contours

    # ptsΣLup = [pts_z1circ[4], z1 + abs(0.5*(z1+pts_z0circ[2]) - z1)*exp(1im*3*pi/4), z0 + abs(0.5*(z1+pts_z0circ[2]) - z1)*exp(1im*phi0left), pts_z0circ[2]]
    # NO LITTLE CIRCLE AT Z1 TEMPORARILY
    ptsΣLup = [Z1, Z1 + abs(0.5*(Z1+ptsZ0circ[2]) - Z1)*exp(1im*3*pi/4), Z0 + abs(0.5*(Z1+ptsZ0circ[2]) - Z1)*exp(1im*ϕ0left), ptsZ0circ[2]]
    ptsΣLdown = conj.(ptsΣLup)|>reverse
    ΣLup = [ptsΣLup[1] ptsΣLup[2]; ptsΣLup[2] ptsΣLup[3]; ptsΣLup[3] ptsΣLup[4]]
    ΣLdown = [ptsΣLdown[1] ptsΣLdown[2]; ptsΣLdown[2] ptsΣLdown[3]; ptsΣLdown[3] ptsΣLdown[4]]
    # The points below may need to start with ptsZ1circ[3] in the first two components.
    ptsΣup = [Z1, Z1+1im*imag(ΣjunctionInsideUHP), ptsZ0circ[3]]
    Σup = [ptsΣup[1] ptsΣup[2]; ptsΣup[2] ptsΣup[3]]
    ptsΣdown = conj.(ptsΣup)|>reverse
    Σdown = [ptsΣdown[1] ptsΣdown[2]; ptsΣdown[2] ptsΣdown[3]]
    # ptsΣRup = [pts_z1circ[2], ΣjunctionPointInsideUHP]
    ptsΣRup = [Z1, ΣjunctionInsideUHP]
    ΣRup = [ptsΣRup[1] ptsΣRup[2]]
    ptsΣRdown = conj.(ptsΣRup)|>reverse
    ΣRdown = [ptsΣRdown[1] ptsΣRdown[2]]
    ptsCommonΣΓtoZ0 = [ptsΣRup[2], ptsZ0circ[4]]
    ptsCommonΣΓfromZ0CC = [ptsZ0CCcirc[2], ptsΣRdown[1]]
    CommonΣΓtoZ0 = [ptsCommonΣΓtoZ0[1] ptsCommonΣΓtoZ0[2]]
    CommonΣΓfromZ0CC = [ptsCommonΣΓfromZ0CC[1] ptsCommonΣΓfromZ0CC[2]]
    
    # I contours

    ptsArcI = [Z1, Z2]
    ArcI = [ptsArcI[1] ptsArcI[2]]


    # Γ contours
    
    ΓjunctionInsideUHP = 0.25*abs(ΣjunctionInsideUHP-Z2)*exp(3im*pi/4) + Z2
    # ΓjunctionInsideLHP = abs(0.5*ΣjunctionInsideUHP+0.5*Z2)*exp(-3im*pi/4) + Z2
    # ptsΓRup may terminate at a little circle point centered at Z2 instead of Z2. 
    ptsΓRup = [ΣjunctionInsideUHP, ΓjunctionInsideUHP, Z2]
    ΓRup = [ptsΓRup[1] ptsΓRup[2]; ptsΓRup[2] ptsΓRup[3]]
    ptsΓRdown = conj.(ptsΓRup)|>reverse
    ΓRdown = [ptsΓRdown[1] ptsΓRdown[2]; ptsΓRdown[2] ptsΓRdown[3]]
    # ptsΓLup may terminate at a little circle point centered at Z2 instead of Z2.
    ptsΓLup = [ptsZ0circ[1], ptsZ0circ[1]+exp(1im*argZ0), Z2+exp(1im*pi/4), Z2]
    ΓLup = [ptsΓLup[1] ptsΓLup[2]; ptsΓLup[2] ptsΓLup[3]; ptsΓLup[3] ptsΓLup[4]]
    ptsΓLdown = conj.(ptsΓLup)|>reverse
    ΓLdown = [ptsΓLdown[1] ptsΓLdown[2]; ptsΓLdown[2] ptsΓLdown[3]; ptsΓLdown[3] ptsΓLdown[4]]

    # Genuz-0 g-function
     # R(z) is cut on Σ and satisfies R(z)^2 = (z-z0)(z-z0cc), r(z)-z = o(1) as z→∞.
    # rbranch is cut on line segments [Z0cc,Z1] and [Z1,Z0]
    Rbranch = z -> sqrt((z-Z0)/(z-Z1)|>ComplexF64)*sqrt((z-Z1)/(z-Z0CC)|>ComplexF64)*(z-Z0CC)
    # rright_vercut the right BV for R(z) cut on line segments [Z0cc,-i∞] and [Z0,+i∞]
    # Rright_vercut = z -> sqrt((z-Z0)*(z-Z0CC)|>ComplexF64)
    # rleft_vercut the left BV for R(z) cut on line segments [Z0cc,-i∞] and [Z0,+i∞]
    # Rleft_vercut = z -> -sqrt((z-Z0)*(z-Z0CC)|>ComplexF64)
    # rright_horcutneginfty is the right BV for R(z) with horizontal cuts to -∞.
    Rright_horcutneginfty = z -> sqrt((z-Z0|>ComplexF64))*sqrt((z-Z0CC|>ComplexF64))
    # rleftt_horcutneginfty is the left BV for r(z) with horizontal cuts to +∞.
    Rleft_horcutposinfty = z -> -sqrt((Z0-z|>ComplexF64))*sqrt((Z0CC-z|>ComplexF64))
    # Rstraight has vertical cut between z0 and z0cc, roundoff is different.
    Rstraight_uhp = z -> ((z-Z0)/(z-Z0CC))^(1/2)*(z-Z0CC)
    Rstraight_lhp = z -> ((z-Z0CC)/(z-Z0))^(1/2)*(z-Z0)

    # The g-function and the h-function cut on Σ.
    # The implementations cut on Σ
    # gbranch = z -> (Rbranch(z))^3/z - θ(z) - 3/2^(1/3) - (w^2)/6
    # hbranch = z -> (Rbranch(z))^3/z - 3/2^(1/3) - (w^2)/6
    # The implementations cut horizontally to ∞ from the branch points
    # gright = z -> (Rright_horcutneginfty(z))^3/z - θ(z) - 3/2^(1/3) - (w^2)/6
    # gleft = z -> (Rleft_horcutposinfty(z)^3)/z - θ(z) - 3/2^(1/3) - (w^2)/6
    
    hright_horcut_neginf = z -> (Rright_horcutneginfty(z)^3)/z - 3.0*2^(-1/3) - (w^2)/6.
    hleft_horcut_posinf = z -> (Rleft_horcutposinfty(z)^3)/z - 3.0*2^(-1/3) - (w^2)/6.
    hstraight_uhp = z -> (Rstraight_uhp(z)^3)/z - 3.0*2^(-1/3) - (w^2)/6.
    hstraight_lhp = z -> (Rstraight_lhp(z)^3)/z - 3.0*2^(-1/3) - (w^2)/6.


    # Jump matrices

    # G = (z->[1 0;0 1])|>mkfun
    # Jumps on Γ
    # Jump matrices relating to ΓLup
    # VL = [1.0 0; -bf/af 1.0]
    # VLi = [1.0 0; bf/af 1.0]
    VLeH = z -> [1.0 0; -(bf/af)*exp(2im*TT*hright_horcut_neginf(z)) 1.0]
    G_ΓLup = mkfun(VLeH)
    J_ΓLup = [G_ΓLup, G_ΓLup, G_ΓLup]
    # Jump matrices relating to ΓRup
    VR = [1.0 af*bf; 0 1.0]
    VRi = [1.0 -af*bf; 0 1.0]
    VReH = z -> [1.0 af*bf*exp(-2im*TT*hright_horcut_neginf(z)); 0 1.0]
    G_ΓRup = mkfun(VReH)
    J_ΓRup = [G_ΓRup, G_ΓRup]
    # Jumps relating to ΓRdown
    YR = [1.0 0; -af*bf 1.0]
    YRi = [1.0 0; af*bf 1.0]
    YReH = z -> [1.0 0; -af*bf*exp(2im*TT*hright_horcut_neginf(z)) 1.0]
    G_ΓRdown = mkfun(YReH)
    J_ΓRdown = [G_ΓRdown, G_ΓRdown]
    # Jumps relating to ΓLdown
    YL = [1.0 bf/af; 0 1.0]
    YLi = [1.0 -bf/af; 0 1.0]
    YLeH = z -> [1.0 (bf/af)*exp(-2im*TT*hright_horcut_neginf(z)); 0 1.0]
    G_ΓLdown = mkfun(YLeH)
    J_ΓLdown = [G_ΓLdown, G_ΓLdown, G_ΓLdown]
    
    # Jump on I
    Ds3 = [af^(-2) 0; 0 af^2]
    Ds3eH = z -> [af^(-2) 0; 0 af^2]
    G_I = mkfun(Ds3eH)
    J_I = [G_I]

    # Jumps relating to ΣLup
    # WL = [1.0 -af/bf; 0 1.0]
    WLi = [1.0 af/bf; 0 1.0]
    WLeH = z -> [1.0 -(af/bf)*exp(-2im*TT*hleft_horcut_posinf(z)); 0 1.0]
    G_ΣLup = mkfun(WLeH)
    J_ΣLup = [G_ΣLup, G_ΣLup, G_ΣLup]
    # Jumps relating to Σup
    W = [0. af/bf; -bf/af 0.]
    Wi = [0. -af/bf; bf/af 0.]
    WeH = z -> [0. (af/bf)*exp(-1im*TT*κ); -(bf/af)*exp(1im*TT*κ) 0.]
    G_Σup = mkfun(WeH)
    J_Σup = [G_Σup, G_Σup]
    # Jumps relating to ΣRup
    # WR = [1.0 -af^3/bf; 0 1.0]
    # WRi = [1.0 af^3/bf; 0 1.0]
    WReH = z -> [1.0 -(af^3/bf)*exp(-2im*TT*hright_horcut_neginf(z)); 0 1.0]
    G_ΣRup = mkfun(WReH)
    J_ΣRup = [G_ΣRup]
    # Jumps relating to ΣLdown
    # XL = [1.0 0; af/bf 1.0]
    XLi = [1.0 0; -af/bf 1.0]
    XLeH = z -> [1.0 0; (af/bf)*exp(2im*TT*hleft_horcut_posinf(z)) 1.0]
    G_ΣLdown = mkfun(XLeH)
    J_ΣLdown = [G_ΣLdown, G_ΣLdown, G_ΣLdown]
    # Jumps relating to Σdown
    # X = [0. bf/af; -af/bf 0.]
    Xi = [0. -bf/af; af/bf 0.]
    XeH = z -> [0. (bf/af)*exp(-1im*TT*κ); -(af/bf)*exp(1im*TT*κ) 0.]
    G_Σdown = mkfun(XeH)
    J_Σdown = [G_Σdown, G_Σdown]
    # Jumps relating to ΣRdown
    # XR = [1.0 0; af^3/bf 1.0]
    # XRi = [1.0 0; -af^3/bf 1.0]
    XReH = z -> [1.0 0; (af^3/bf)*exp(2im*TT*hright_horcut_neginf(z)) 1.0]
    G_ΣRdown = mkfun(XReH)
    J_ΣRdown = [G_ΣRdown]
    # Collapsed product jumps Pup and Pdown.
    # Pup = VRi*WR
    Pup = [1.0 -af/bf; 0 1.0]
    PupeH = z -> [1.0 -(af/bf)*exp(-2im*TT*hright_horcut_neginf(z)); 0 1.0]
    J_commonup = [mkfun(PupeH)]
    # Pdown = YRi*XR
    Pdown = [1.0 0; af/bf 1.0]
    PdowneH = z -> [1.0 0; (af/bf)*exp(2im*TT*hright_horcut_neginf(z)) 1.0]
    J_commondown = [mkfun(PdowneH)]
    # return(Pdown-(YRi*XR), Pup-(VRi*WR))

    # Jump matrices relating to Z0circ
    G_Z0_12 = (z->[exp(-1im*TT*hstraight_uhp(z)) 0; 0 exp(1im*TT*hstraight_uhp(z))]*WLi*Wi)|>mkfun
    G_Z0_23 = (z->[exp(-1im*TT*hleft_horcut_posinf(z)) 0; 0 exp(1im*TT*hleft_horcut_posinf(z))]*Wi)|>mkfun
    G_Z0_34 = (z->[exp(-1im*TT*hright_horcut_neginf(z)) 0; 0 exp(1im*TT*hright_horcut_neginf(z))])|>mkfun
    G_Z0_41 = (z->[exp(-1im*TT*hright_horcut_neginf(z)) 0; 0 exp(1im*TT*hright_horcut_neginf(z))]*Pup)|>mkfun
    J_Z0 = [G_Z0_12, G_Z0_23, G_Z0_34, G_Z0_41]

    # Jump matrices relating to Z0CCcirc
    G_Z0CC_12 = (z->[exp(-1im*TT*hright_horcut_neginf(z)) 0; 0 exp(1im*TT*hright_horcut_neginf(z))]*Pdown)|>mkfun
    G_Z0CC_23 = (z->[exp(-1im*TT*hright_horcut_neginf(z)) 0; 0 exp(1im*TT*hright_horcut_neginf(z))])|>mkfun
    G_Z0CC_34 = (z->[exp(-1im*TT*hleft_horcut_posinf(z)) 0; 0 exp(1im*TT*hleft_horcut_posinf(z))]*Xi)|>mkfun
    G_Z0CC_41 = (z->[exp(-1im*TT*hleft_horcut_posinf(z)) 0; 0 exp(1im*TT*hleft_horcut_posinf(z))]*XLi*Xi)|>mkfun
    J_Z0CC = [G_Z0CC_12, G_Z0CC_23, G_Z0CC_34, G_Z0CC_41]
    # J_Z0 = [G, G, G, G]
    # J_Z0CC = [G, G, G, G]
    # J_ΣLup = [G, G, G]
    # J_Σup = [G, G]
    # J_ΣRup = [G]
    # J_commonup = [G]
    # J_ΓRup = [G, G]
    # J_ΓLup = [G, G, G]

    # J_ΣLdown = [G, G, G]
    # J_Σdown = [G, G]
    # J_ΣRdown = [G]
    # J_commondown = [G]
    # J_ΓRdown = [G, G]
    # J_ΓLdown = [G, G, G]

    # Put together the arcs and jump matrices. Construct the RHP.
    circles = vcat(Z0circ, Z0CCcirc)
    ΣΓup = vcat(ΣLup, Σup, ΣRup, CommonΣΓtoZ0, ΓRup, ΓLup)
    ΣΓdown = vcat(ΣLdown, Σdown, ΣRdown, CommonΣΓfromZ0CC, ΓRdown, ΓLdown)
    intervals = vcat(circles, ΣΓup, ΣΓdown, ArcI)

    J = vcat(J_Z0, J_Z0CC, J_ΣLup, J_Σup, J_ΣRup, J_commonup, J_ΓRup, J_ΓLup, J_ΣLdown, J_Σdown, J_ΣRdown, J_commondown, J_ΓRdown, J_ΓLdown, J_I)

    rhp = RHP(intervals, J)
    rhp_solver = RHSolver(rhp)
    numarcs = length(J)
    u1 = rhp_solver([1 0],numarcs*npts)
    u2 = rhp_solver([0 1],numarcs*npts)
    U = rhp_solver(([1 0], [0 1]), numarcs*npts)
    U = hcat(U...) |> transpose |> Matrix
    Usums = sum.(U)
    nls = -(1/π)*Usums[1,2]
    nls = nls*exp(-1im*angle(a*b))/TT
    return nls, rhp, rhp_solver, U, u1, u2
end


function rwio_largeT(T,w,a,b, npts::Integer)
    TT = T^(1/3)
    # θ = z -> w*z+z^2+2/z
    if w >= WCRIT
        println("w value is too large.")
        return
    end
    af = abs(a)/sqrt(abs(a)^2 + abs(b)^2)
    bf = abs(b)/sqrt(abs(a)^2 + abs(b)^2)
    # τ = bf/af
    Z1 = (1/2)*(-w/6 - sqrt(w^2/36+2^(5/3)))
    Z2 = (1/2)*(-w/6 + sqrt(w^2/36+2^(5/3)))
    Z0 = (1/3)*(-w + 1im*sqrt(WCRIT^2-w^2))
    Z0CC = (1/3)*(-w - 1im*sqrt(WCRIT^2-w^2))
    # argZ0Z1 = angle(Z0-Z1)
    # midZ1Z2 = 0.5*(Z1+Z2)
    # κ is the constant sum of the boundary values of h(z) on Σ.
    κ = -108^(1/3) - (w^2)/3
    # Rotating angle for the polygon at z0.
    argZ0 = angle_z0z1(w)
     # Radii of small disks
    # radZ0 = min(0.5/(1+abs(T)^(2/9)), 0.3*abs(Z1-Z0), 0.5)
    radZ0 = min(2.5/(1+abs(T)^(2/9)), 0.3*abs(Z1-Z0), 0.25)
    # radZ1 = min(2/(1+abs(T)^(1/6)), 0.3*abs(Z1-Z0), 0.25)
    # radZ1 = min(0.5/(1+abs(T)^(1/6)), 0.3*abs(Z1-Z0), 0.5)

    # Small circle points
    ptsZ0circ = [1,exp(2im*pi/3),exp(1im*pi),exp(-2im*pi/3)]*radZ0*exp(1im*argZ0) .+ Z0
    Z0circ = [ptsZ0circ[1] ptsZ0circ[2]; ptsZ0circ[2] ptsZ0circ[3]; ptsZ0circ[3] ptsZ0circ[4]; ptsZ0circ[4] ptsZ0circ[1]]
    ptsZ0CCcirc = [1,exp(2im*pi/3),exp(1im*pi),exp(-2im*pi/3)]*radZ0*exp(-1im*argZ0) .+ Z0CC
    Z0CCcirc = [ptsZ0CCcirc[1] ptsZ0CCcirc[2]; ptsZ0CCcirc[2] ptsZ0CCcirc[3]; ptsZ0CCcirc[3] ptsZ0CCcirc[4]; ptsZ0CCcirc[4] ptsZ0CCcirc[1]]
    # ptsZ1circ = [1, exp(1im*pi/4), exp(1im*pi/2), exp(3im*pi/4), exp(-3im*pi/4), exp(-1im*pi/2), exp(-1im*pi/4), 1]*radZ1 .+ Z1
    # ptsZ2circ = [exp(1im*pi), exp(-3im*pi/4), exp(-1im*pi/4), exp(1im*pi/4), exp(3im*pi/4)]*radZ1 .+ Z2

    # Internal parameters required for the contours Σ and Γ
    # Exit angle from z0 to the exterior domain
    ϕ0left = angle(ptsZ0circ[2]-Z0)
    # Exit angle from z0 to the interior domain
    ϕ0right = angle(ptsZ0circ[4]-Z0)
    # t1 is the parameter so that Z1+t1*exp(1im*pi/4) intersects the line coming out of Z0 to the interior domain with the Airy angle.
    t1 = (-1/sin(ϕ0right-pi/4))*(imag(Z0)*cos(ϕ0right) + (Z1-real(Z0))*sin(ϕ0right))
    # The midwert below may start with ptsZ1circ instead of Z1 in case there is a little circle at Z1
    # midvert = imag(0.5*(ptsZ1circ[2]+ptsZ0circ[3]))
    # midvert = imag(0.5*(Z1+ptsZ0circ[3]))
    # ΣjunctionInsideUHP = midZ1Z2+1im*midvert
    # ΣjunctionInsideLHP = midZ1Z2-1im*midvert
    ΣjunctionInsideUHP = Z1 + t1*exp(1im*pi/4)
    # ΣjunctionInsideLHP = Z1 + t1*exp(-1im*pi/4)

    # Σ contours

    # ptsΣLup = [pts_z1circ[4], z1 + abs(0.5*(z1+pts_z0circ[2]) - z1)*exp(1im*3*pi/4), z0 + abs(0.5*(z1+pts_z0circ[2]) - z1)*exp(1im*phi0left), pts_z0circ[2]]
    # NO LITTLE CIRCLE AT Z1 TEMPORARILY
    ptsΣLup = [Z1, Z1 + abs(0.5*(Z1+ptsZ0circ[2]) - Z1)*exp(1im*3*pi/4), Z0 + abs(0.5*(Z1+ptsZ0circ[2]) - Z1)*exp(1im*ϕ0left), ptsZ0circ[2]]
    ptsΣLdown = conj.(ptsΣLup)|>reverse
    ΣLup = [ptsΣLup[1] ptsΣLup[2]; ptsΣLup[2] ptsΣLup[3]; ptsΣLup[3] ptsΣLup[4]]
    ΣLdown = [ptsΣLdown[1] ptsΣLdown[2]; ptsΣLdown[2] ptsΣLdown[3]; ptsΣLdown[3] ptsΣLdown[4]]
    # The points below may need to start with ptsZ1circ[3] in the first two components.
    ptsΣup = [Z1, Z1+1im*imag(ΣjunctionInsideUHP), ptsZ0circ[3]]
    Σup = [ptsΣup[1] ptsΣup[2]; ptsΣup[2] ptsΣup[3]]
    ptsΣdown = conj.(ptsΣup)|>reverse
    Σdown = [ptsΣdown[1] ptsΣdown[2]; ptsΣdown[2] ptsΣdown[3]]
    # ptsΣRup = [pts_z1circ[2], ΣjunctionPointInsideUHP]
    ptsΣRup = [Z1, ΣjunctionInsideUHP]
    ΣRup = [ptsΣRup[1] ptsΣRup[2]]
    ptsΣRdown = conj.(ptsΣRup)|>reverse
    ΣRdown = [ptsΣRdown[1] ptsΣRdown[2]]
    ptsCommonΣΓtoZ0 = [ptsΣRup[2], ptsZ0circ[4]]
    ptsCommonΣΓfromZ0CC = [ptsZ0CCcirc[2], ptsΣRdown[1]]
    CommonΣΓtoZ0 = [ptsCommonΣΓtoZ0[1] ptsCommonΣΓtoZ0[2]]
    CommonΣΓfromZ0CC = [ptsCommonΣΓfromZ0CC[1] ptsCommonΣΓfromZ0CC[2]]
    
    # I contours

    ptsArcI = [Z1, Z2]
    ArcI = [ptsArcI[1] ptsArcI[2]]


    # Γ contours
    
    ΓjunctionInsideUHP = 0.25*abs(ΣjunctionInsideUHP-Z2)*exp(3im*pi/4) + Z2
    # ΓjunctionInsideLHP = abs(0.4*ΣjunctionInsideUHP+0.6*Z2)*exp(-3im*pi/4) + Z2
    # ptsΓRup may terminate at a little circle point centered at Z2 instead of Z2. 
    ptsΓRup = [ΣjunctionInsideUHP, ΓjunctionInsideUHP, Z2]
    ΓRup = [ptsΓRup[1] ptsΓRup[2]; ptsΓRup[2] ptsΓRup[3]]
    ptsΓRdown = conj.(ptsΓRup)|>reverse
    ΓRdown = [ptsΓRdown[1] ptsΓRdown[2]; ptsΓRdown[2] ptsΓRdown[3]]
    # ptsΓLup may terminate at a little circle point centered at Z2 instead of Z2.
    ptsΓLup = [ptsZ0circ[1], ptsZ0circ[1]+exp(1im*argZ0), Z2+exp(1im*pi/4), Z2]
    ΓLup = [ptsΓLup[1] ptsΓLup[2]; ptsΓLup[2] ptsΓLup[3]; ptsΓLup[3] ptsΓLup[4]]
    ptsΓLdown = conj.(ptsΓLup)|>reverse
    ΓLdown = [ptsΓLdown[1] ptsΓLdown[2]; ptsΓLdown[2] ptsΓLdown[3]; ptsΓLdown[3] ptsΓLdown[4]]

    # Genuz-0 g-function
     # R(z) is cut on Σ and satisfies R(z)^2 = (z-z0)(z-z0cc), r(z)-z = o(1) as z→∞.
    # rbranch is cut on line segments [Z0cc,Z1] and [Z1,Z0]
    Rbranch = z -> sqrt((z-Z0)/(z-Z1)|>ComplexF64)*sqrt((z-Z1)/(z-Z0CC)|>ComplexF64)*(z-Z0CC)
    # rright_vercut the right BV for R(z) cut on line segments [Z0cc,-i∞] and [Z0,+i∞]
    # Rright_vercut = z -> sqrt((z-Z0)*(z-Z0CC)|>ComplexF64)
    # rleft_vercut the left BV for R(z) cut on line segments [Z0cc,-i∞] and [Z0,+i∞]
    # Rleft_vercut = z -> -sqrt((z-Z0)*(z-Z0CC)|>ComplexF64)
    # rright_horcutneginfty is the right BV for R(z) with horizontal cuts to -∞.
    Rright_horcutneginfty = z -> sqrt((z-Z0|>ComplexF64))*sqrt((z-Z0CC|>ComplexF64))
    # rleftt_horcutneginfty is the left BV for r(z) with horizontal cuts to +∞.
    Rleft_horcutposinfty = z -> -sqrt((Z0-z|>ComplexF64))*sqrt((Z0CC-z|>ComplexF64))
    # Rstraight has vertical cut between z0 and z0cc, roundoff is different.
    Rstraight_uhp = z -> ((z-Z0)/(z-Z0CC))^(1/2)*(z-Z0CC)
    Rstraight_lhp = z -> ((z-Z0CC)/(z-Z0))^(1/2)*(z-Z0)

    # The g-function and the h-function cut on Σ.
    # The implementations cut on Σ
    # gbranch = z -> (Rbranch(z))^3/z - θ(z) - 3/2^(1/3) - (w^2)/6
    # hbranch = z -> (Rbranch(z))^3/z - 3/2^(1/3) - (w^2)/6
    # The implementations cut horizontally to ∞.
    # gright = z -> (Rright_horcutneginfty(z))^3/z - θ(z) - 3/2^(1/3) - (w^2)/6
    # gleft = z -> (Rleft_horcutposinfty(z)^3)/z - θ(z) - 3/2^(1/3) - (w^2)/6
    
    hright_horcut_neginf = z -> (Rright_horcutneginfty(z)^3)/z - 3.0*2^(-1/3) - (w^2)/6.
    hleft_horcut_posinf = z -> (Rleft_horcutposinfty(z)^3)/z - 3.0*2^(-1/3) - (w^2)/6.
    hstraight_uhp = z -> (Rstraight_uhp(z)^3)/z - 3.0*2^(-1/3) - (w^2)/6.
    hstraight_lhp = z -> (Rstraight_lhp(z)^3)/z - 3.0*2^(-1/3) - (w^2)/6.


    # Jump matrices

    # G = (z->[1 0;0 1])|>mkfun
    # Jumps on Γ
    # Jump matrices relating to ΓLup
    # VL = [1.0 0; -bf/af 1.0]
    # VLi = [1.0 0; bf/af 1.0]
    VLeH = z -> [1.0 0; -(bf/af)*exp(2im*TT*hright_horcut_neginf(z)) 1.0]
    G_ΓLup = mkfun(VLeH)
    J_ΓLup = [G_ΓLup, G_ΓLup, G_ΓLup]
    # Jump matrices relating to ΓRup
    # VR = [1.0 af*bf; 0 1.0]
    # VRi = [1.0 -af*bf; 0 1.0]
    VReH = z -> [1.0 af*bf*exp(-2im*TT*hright_horcut_neginf(z)); 0 1.0]
    G_ΓRup = mkfun(VReH)
    J_ΓRup = [G_ΓRup, G_ΓRup]
    # Jumps relating to ΓRdown
    # YR = [1.0 0; -af*bf 1.0]
    # YRi = [1.0 0; af*bf 1.0]
    YReH = z -> [1.0 0; -af*bf*exp(2im*TT*hright_horcut_neginf(z)) 1.0]
    G_ΓRdown = mkfun(YReH)
    J_ΓRdown = [G_ΓRdown, G_ΓRdown]
    # Jumps relating to ΓLdown
    # YL = [1.0 bf/af; 0 1.0]
    # YLi = [1.0 -bf/af; 0 1.0]
    YLeH = z -> [1.0 (bf/af)*exp(-2im*TT*hright_horcut_neginf(z)); 0 1.0]
    G_ΓLdown = mkfun(YLeH)
    J_ΓLdown = [G_ΓLdown, G_ΓLdown, G_ΓLdown]
    
    # Jump on I
    # Ds3 = [af^(-2) 0; 0 af^2]
    Ds3eH = z -> [af^(-2) 0; 0 af^2]
    G_I = mkfun(Ds3eH)
    J_I = [G_I]

    # Jumps relating to ΣLup
    # WL = [1.0 -af/bf; 0 1.0]
    WLi = [1.0 af/bf; 0 1.0]
    WLeH = z -> [1.0 -(af/bf)*exp(-2im*TT*hleft_horcut_posinf(z)); 0 1.0]
    G_ΣLup = mkfun(WLeH)
    J_ΣLup = [G_ΣLup, G_ΣLup, G_ΣLup]
    # Jumps relating to Σup
    # W = [0. af/bf; -bf/af 0.]
    Wi = [0. -af/bf; bf/af 0.]
    WeH = z -> [0. (af/bf)*exp(-1im*TT*κ); -(bf/af)*exp(1im*TT*κ) 0.]
    G_Σup = mkfun(WeH)
    J_Σup = [G_Σup, G_Σup]
    # Jumps relating to ΣRup
    # WR = [1.0 -af^3/bf; 0 1.0]
    # WRi = [1.0 af^3/bf; 0 1.0]
    WReH = z -> [1.0 -(af^3/bf)*exp(-2im*TT*hright_horcut_neginf(z)); 0 1.0]
    G_ΣRup = mkfun(WReH)
    J_ΣRup = [G_ΣRup]
    # Jumps relating to ΣLdown
    # XL = [1.0 0; af/bf 1.0]
    XLi = [1.0 0; -af/bf 1.0]
    XLeH = z -> [1.0 0; (af/bf)*exp(2im*TT*hleft_horcut_posinf(z)) 1.0]
    G_ΣLdown = mkfun(XLeH)
    J_ΣLdown = [G_ΣLdown, G_ΣLdown, G_ΣLdown]
    # Jumps relating to Σdown
    # X = [0. bf/af; -af/bf 0.]
    Xi = [0. -bf/af; af/bf 0.]
    XeH = z -> [0. (bf/af)*exp(-1im*TT*κ); -(af/bf)*exp(1im*TT*κ) 0.]
    G_Σdown = mkfun(XeH)
    J_Σdown = [G_Σdown, G_Σdown]
    # Jumps relating to ΣRdown
    # XR = [1.0 0; af^3/bf 1.0]
    # XRi = [1.0 0; -af^3/bf 1.0]
    XReH = z -> [1.0 0; (af^3/bf)*exp(2im*TT*hright_horcut_neginf(z)) 1.0]
    G_ΣRdown = mkfun(XReH)
    J_ΣRdown = [G_ΣRdown]
    # Collapsed product jumps Pup and Pdown.
    # Pup = VRi*WR
    Pup = [1.0 -af/bf; 0 1.0]
    PupeH = z -> [1.0 -(af/bf)*exp(-2im*TT*hright_horcut_neginf(z)); 0 1.0]
    J_commonup = [mkfun(PupeH)]
    # Pdown = YRi*XR
    Pdown = [1.0 0; af/bf 1.0]
    PdowneH = z -> [1.0 0; (af/bf)*exp(2im*TT*hright_horcut_neginf(z)) 1.0]
    J_commondown = [mkfun(PdowneH)]
    # return(Pdown-(YRi*XR), Pup-(VRi*WR))

    # Jump matrices relating to Z0circ
    G_Z0_12 = (z->[exp(-1im*TT*hstraight_uhp(z)) 0; 0 exp(1im*TT*hstraight_uhp(z))]*WLi*Wi)|>mkfun
    G_Z0_23 = (z->[exp(-1im*TT*hleft_horcut_posinf(z)) 0; 0 exp(1im*TT*hleft_horcut_posinf(z))]*Wi)|>mkfun
    G_Z0_34 = (z->[exp(-1im*TT*hright_horcut_neginf(z)) 0; 0 exp(1im*TT*hright_horcut_neginf(z))])|>mkfun
    G_Z0_41 = (z->[exp(-1im*TT*hright_horcut_neginf(z)) 0; 0 exp(1im*TT*hright_horcut_neginf(z))]*Pup)|>mkfun
    J_Z0 = [G_Z0_12, G_Z0_23, G_Z0_34, G_Z0_41]

    # Jump matrices relating to Z0CCcirc
    G_Z0CC_12 = (z->[exp(-1im*TT*hright_horcut_neginf(z)) 0; 0 exp(1im*TT*hright_horcut_neginf(z))]*Pdown)|>mkfun
    G_Z0CC_23 = (z->[exp(-1im*TT*hright_horcut_neginf(z)) 0; 0 exp(1im*TT*hright_horcut_neginf(z))])|>mkfun
    G_Z0CC_34 = (z->[exp(-1im*TT*hleft_horcut_posinf(z)) 0; 0 exp(1im*TT*hleft_horcut_posinf(z))]*Xi)|>mkfun
    G_Z0CC_41 = (z->[exp(-1im*TT*hleft_horcut_posinf(z)) 0; 0 exp(1im*TT*hleft_horcut_posinf(z))]*XLi*Xi)|>mkfun
    J_Z0CC = [G_Z0CC_12, G_Z0CC_23, G_Z0CC_34, G_Z0CC_41]
    # J_Z0 = [G, G, G, G]
    # J_Z0CC = [G, G, G, G]
    # J_ΣLup = [G, G, G]
    # J_Σup = [G, G]
    # J_ΣRup = [G]
    # J_commonup = [G]
    # J_ΓRup = [G, G]
    # J_ΓLup = [G, G, G]

    # J_ΣLdown = [G, G, G]
    # J_Σdown = [G, G]
    # J_ΣRdown = [G]
    # J_commondown = [G]
    # J_ΓRdown = [G, G]
    # J_ΓLdown = [G, G, G]

    # Put together the arcs and jump matrices. Construct the RHP.
    circles = vcat(Z0circ, Z0CCcirc)
    ΣΓup = vcat(ΣLup, Σup, ΣRup, CommonΣΓtoZ0, ΓRup, ΓLup)
    ΣΓdown = vcat(ΣLdown, Σdown, ΣRdown, CommonΣΓfromZ0CC, ΓRdown, ΓLdown)
    intervals = vcat(circles, ΣΓup, ΣΓdown, ArcI)

    J = vcat(J_Z0, J_Z0CC, J_ΣLup, J_Σup, J_ΣRup, J_commonup, J_ΓRup, J_ΓLup, J_ΣLdown, J_Σdown, J_ΣRdown, J_commondown, J_ΓRdown, J_ΓLdown, J_I)

    rhp = RHP(intervals, J)
    rhp_solver = RHSolver(rhp)
    numarcs = length(J)
    # u1 = rhp_solver([1 0],numarcs*npts)
    # u2 = rhp_solver([0 1],numarcs*npts)
    U = rhp_solver(([1 0], [0 1]), numarcs*npts)
    GC.gc()
    U = hcat(U...) |> transpose |> Matrix
    Usums = sum.(U)
    nls = -(1/π)*Usums[1,2]
    nls = nls*exp(-1im*angle(a*b))/TT
    return nls
end

# end # module

