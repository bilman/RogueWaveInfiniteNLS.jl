# SPDX-License-Identifier: MIT
module RogueWaveInfiniteNLS
using LinearAlgebra
using OperatorApproximation
using Polynomials
export AbstractPolynomial, AbstractUnivariatePolynomial, ArnoldiFit, ChebyshevT, FactoredPolynomial, ImmutablePolynomial, LaurentPolynomial, Polynomial, Polynomials, RationalFunction, SparsePolynomial, chop!, coeffs, companion, degree, derivative, fit, fromroots, hasnan, integrate, isintegral, ismonic, lowest_terms, mapdomain, poles, printpoly, residues, roots, truncate!, vander, variable

# Export basic user-end functions
include("rwio_RHP.jl")
include("rwio.jl")



greet() = print("Hello World!")

end # module RogueWaveInfiniteNLS
