# SPDX-License-Identifier: MIT
module RogueWaveInfiniteNLS

using LinearAlgebra
using OperatorApproximation
using Polynomials

# Export basic user-end functions
include("rwio_RHP.jl")
include("rwio.jl")

export psi, psi_undeformed, psi_largeX, psi_largeT, psi_Painleve, Psi_slice_T0
export rwio_undeformed, rwio_largeX, rwio_largeT, rwio_Painleve, rwio_undeformed_RHP, rwio_largeX_RHP, rwio_largeT_RHP, rwio_Painleve_RHP
export vfromTw, vfromXT, wfromXT, Tmax_givenX, Xmax_givenT, XfromTw, TfromXv
export NODEF_PTS, LARGEX_PTS, LARGET_PTS, VCRIT, WCRIT

end # module RogueWaveInfiniteNLS
