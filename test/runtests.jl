using RogueWaveInfiniteNLS
using Test

# Check if the solution at the peak is computed correctly
@testset "RogueWaveInfiniteNLS.jl: Peak test" begin
    atest = 1
    btest = 1
    truth = 8*conj(atest)*conj(btest)/(abs(atest)^2+abs(btest)^2)
    @test abs(psi(0.,0., atest, btest, 1) - truth) < 1e-10

    atest = 1
    btest = 2
    truth = 8*conj(atest)*conj(btest)/(abs(atest)^2+abs(btest)^2)
    @test abs(psi(0.,0., atest, btest, 1) - truth) < 1e-10
    
    atest = 1
    btest = 2+1im
    truth = 8*conj(atest)*conj(btest)/(abs(atest)^2+abs(btest)^2)
    @test abs(psi(0.,0., atest, btest, 1) - truth) < 1e-10

    atest = -1
    btest = 1
    truth = 8*conj(atest)*conj(btest)/(abs(atest)^2+abs(btest)^2)
    @test abs(psi(0.,0., atest, btest, 1) - truth) < 1e-10
end

@testset "RogueWaveInfiniteNLS.jl: Large-X deformations test" begin
    atest = 1.0
    btest = 1.0
    Xtest = 2.5
    Ttest = 0.15
    vtest = vfromXT(Xtest,Ttest)
    truth = rwio_undeformed(Xtest,Ttest, atest, btest, NODEF_PTS)
    @test abs(rwio_largeX(Xtest, vtest, atest, btest, LARGEX_PTS) - truth) < 1e-10

    atest = 1.0
    btest = -1.0+3.0im
    Xtest = 2.5
    Ttest = 0.15
    vtest = vfromXT(Xtest,Ttest)
    truth = rwio_undeformed(Xtest,Ttest, atest, btest, NODEF_PTS)
    @test abs(rwio_largeX(Xtest, vtest, atest, btest, LARGEX_PTS) - truth) < 1e-10
end


@testset "RogueWaveInfiniteNLS.jl: Large-T deformations test" begin
    atest = 1.0
    btest = 1.0
    Xtest = 0.5
    Ttest = 1.2
    wtest = wfromXT(Xtest,Ttest)
    truth = rwio_undeformed(Xtest,Ttest, atest, btest, NODEF_PTS)
    @test abs(rwio_largeT(Ttest, wtest, atest, btest, LARGET_PTS) - truth) < 1e-10

    atest = 1.0
    btest = -1.0+3.0im
    Xtest = 0.5
    Ttest = 1.2
    vtest = vfromXT(Xtest,Ttest)
    truth = rwio_undeformed(Xtest,Ttest, atest, btest, NODEF_PTS)
    @test abs(rwio_largeT(Ttest, wtest, atest, btest, LARGET_PTS) - truth) < 1e-10
end

@testset "RogueWaveInfiniteNLS.jl: Painlevé deformations test" begin
    atest = 1.0
    btest = 1.0
    vtest = VCRIT
    Ttest = 1.2
    Xtest = XfromTw(Ttest, WCRIT)
    truth = rwio_undeformed(Xtest,Ttest, atest, btest, NODEF_PTS)
    @test abs(rwio_Painleve(Xtest, vtest, atest, btest, LARGEX_PTS) - truth) < 1e-10

    atest = 1.0
    btest = -1.0+3.0im
    vtest = VCRIT
    Ttest = 1.2
    Xtest = XfromTw(Ttest, WCRIT)
    truth = rwio_undeformed(Xtest,Ttest, atest, btest, NODEF_PTS)
    @test abs(rwio_Painleve(Xtest, vtest, atest, btest, LARGEX_PTS) - truth) < 1e-10
end

@testset "RogueWaveInfiniteNLS.jl: Painlevé vs Large-X test" begin
    atest = 1.0
    btest = -1.0+3.0im
    vtest = VCRIT*0.99
    Xtest = 400
    Ttest = TfromXv(Xtest, vtest)
    psival_Pain = rwio_Painleve(Xtest,vtest, atest, btest, LARGEX_PTS)
    psival_largeX = rwio_largeX(Xtest, vtest, atest, btest, LARGEX_PTS)
    @test abs(psival_Pain - psival_largeX) < 1e-10
end

@testset "RogueWaveInfiniteNLS.jl: Painlevé vs Large-T test" begin
    atest = 1.0-2.0im
    btest = 3.0+4.0im
    Ttest = 400
    wtest = WCRIT*0.95
    Xtest=XfromTw(Ttest,wtest)
    vtest=vfromXT(Xtest,Ttest)
    psival_Pain = rwio_Painleve(Xtest,vtest, atest, btest, LARGEX_PTS)
    psival_largeT = rwio_largeT(Ttest, wtest, atest, btest, LARGET_PTS)
    @test abs(psival_Pain - psival_largeT) < 1e-10
end

@testset "RogueWaveInfiniteNLS.jl: The symmetry used in the 2nd quadrant" begin
    atest = 2.0+4.0im
    btest = -1.0+3.0im
    Xtest = -1.
    Ttest = 0.35
    truth = rwio_undeformed(Xtest, Ttest, atest, btest, NODEF_PTS)
    testval = psi(Xtest, Ttest, atest, btest, 1)
    @test abs(testval - truth) < 1e-10
end


@testset "RogueWaveInfiniteNLS.jl: The symmetry used in the 3rd quadrant" begin
    atest = 2.0+4.0im
    btest = -1.0+3.0im
    Xtest = -1.
    Ttest = -0.35
    truth = rwio_undeformed(Xtest, Ttest, atest, btest, NODEF_PTS)
    testval = psi(Xtest, Ttest, atest, btest, 1)
    @test abs(testval - truth) < 1e-10
end

@testset "RogueWaveInfiniteNLS.jl: The symmetry used in the 4th quadrant" begin
    atest = 2.0+4.0im
    btest = -1.0+3.0im
    Xtest = 1.
    Ttest = -0.35
    truth = rwio_undeformed(Xtest, Ttest, atest, btest, NODEF_PTS)
    testval = psi(Xtest, Ttest, atest, btest, 1)
    @test abs(testval - truth) < 1e-10
end