# SPDX-License-Identifier: MIT
module RogueWaveInfiniteNLS

using LinearAlgebra
using OperatorApproximation
using Polynomials

# Export basic user-end functions
include("rwio_RHP.jl")
include("rwio.jl")

export psi, psi_nodeformation_rescaled, psi_largeX, psi_largeT, psi_Painleve, rwio_nodeformation_rescaled, rwio_largeX, rwio_largeT, rwio_Painleve, rwio_nodeformation_rescaled_RHP, rwio_largeX_RHP, rwio_largeT_RHP, rwio_Painleve_RHP



# greet() = print("Hello World!")

end # module RogueWaveInfiniteNLS
