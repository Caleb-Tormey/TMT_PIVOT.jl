using TMT_PIVOT
using Test

@testset "TMT_PIVOT.jl" begin
    # Write your tests here.
    @test TMT_PIVOT.add(1, 2) == 3
    @test TMT_PIVOT.greet_TMT_PIVOT() == "Hello TMT_PIVOT!"
    @test TMT_PIVOT.greet_TMT_PIVOT() != "Hello World!"
end
