# SPDX-License-Identifier: MIT
module RogueWaveInfiniteNLS
using LinearAlgebra
using OperatorApproximation
export AbstractZeroOperator, ArgNum, Axis, Basis, BasisExpansion, BlockAbstractOperator, BlockDiagonalAbstractOperator, BoundaryValue, CauchyOperator, CauchyTransform, ChebyshevInterval, ChebyshevMappedInterval, CoefConversion, ConcreteOperator, Conversion, Derivative, DirectSum, DirectedLLobattoMappedInterval, DirectedLobattoMappedInterval, DirectedRLobattoMappedInterval, Domain, Erf, Exterior, FixedGridValues, FloquetDerivative, Fourier, GridAxis, GridDomain, GridValues, Hardy, HermiteFun, HermitePoly, HermiteRealAxis, Interior, Jacobi, JacobiInterval, JacobiMappedInterval, Laurent, Legendre, LobattoInterval, LobattoMappedInterval, MappedCircle, MappedInterval, Multiplication, OperatorApproximation, PeriodicCircle, PeriodicInterval, PeriodicMappedCircle, PeriodicMappedInterval, RHP, RHSolver, RHSolverVec, RealAxis, UltraInterval, UltraMappedInterval, Ultraspherical, UnitCircle, UnitInterval, ZeroOperator, \, adapt, arclength, clearCauchycache, coefplot, coefplot!, dilog, domainplot, domainplot!, eigen, matrix2BlockOperator, mofeval, mult2x2, mvf2mof, plot, plot!, ploteval, ploteval!, rhdomain, rhmult, rhplot, rhrange, rhrhs, rhwellposed, setN, setbasis, ⊕, ⊘, ⊞
using Polynomials
export AbstractPolynomial, AbstractUnivariatePolynomial, ArnoldiFit, ChebyshevT, FactoredPolynomial, ImmutablePolynomial, LaurentPolynomial, Polynomial, Polynomials, RationalFunction, SparsePolynomial, chop!, coeffs, companion, degree, derivative, fit, fromroots, hasnan, integrate, isintegral, ismonic, lowest_terms, mapdomain, poles, printpoly, residues, roots, truncate!, vander, variable

# Export basic user-end functions
include("rwio_RHP.jl")
include("rwio.jl")

export psi, psi_nodeformation_rescaled, psi_largeX, psi_largeT, psi_Painleve, rwio_nodeformation_rescaled, rwio_largeX, rwio_largeT, rwio_Painleve, rwio_nodeformation_rescaled_RHP, rwio_largeX_RHP, rwio_largeT_RHP, rwio_Painleve_RHP



# greet() = print("Hello World!")

end # module RogueWaveInfiniteNLS
