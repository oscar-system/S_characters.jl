using Test
using Oscar

@testset "all S-characters" begin
  # real character table is rational
  tbl = character_table("S4")
  res = s_characters(tbl, rational = false, ppow_nonzero = false)
  @test length(res) == 24
  @test res == s_characters(tbl, rational = false, ppow_nonzero = false, irrats = :qqbar)
  @test res == s_characters(tbl, rational = true, ppow_nonzero = false)

  # real character table is not rational
  tbl = character_table("A5")
  res = s_characters(tbl, rational = false, ppow_nonzero = false)
  @test length(res) == 24
  @test res == s_characters(tbl, rational = false, ppow_nonzero = false, irrats = :qqbar)
end

@testset "only S-characters that are nonzero on all p-power elements" begin
  # real character table is rational
  tbl = character_table("A8")
  res = s_characters(tbl, rational = false, ppow_nonzero = true)
  @test length(res) == 1
  @test res == s_characters(tbl, rational = false, ppow_nonzero = true, irrats = :qqbar)
  @test res == s_characters(tbl, rational = true, ppow_nonzero = true, irrats = :qqbar)

  # real character table is not rational
  tbl = character_table("J1")
  res = s_characters(tbl, rational = false, ppow_nonzero = true)
  @test length(res) == 1
# The computations with `irrats = :qqbar` need quite some time.
# @test res == s_characters(tbl, rational = false, ppow_nonzero = true, irrats = :qqbar)
  @test res == s_characters(tbl, rational = true, ppow_nonzero = true, irrats = :qqbar)
end

@testset "two-step computation" begin
  tbl = character_table("2.A8")
  res = s_characters(tbl, rational = false, ppow_nonzero = true)
  @test length(res) == 3
  @test res == s_characters(tbl, rational = true, ppow_nonzero = true)
  @test sort(res, by = degree) == sort(s_characters_via_factor_group(tbl, [1, 2], rational = true, ppow_nonzero = true), by = degree)
  @test sort(res, by = degree) == sort(s_characters_via_factor_group(tbl, [1, 2], rational = false, ppow_nonzero = true), by = degree)
end

@testset "the tests from the docstrings" begin
  tbl = character_table("A5");
  P, galoissums, ppow_pos = s_character_simplex(tbl);
  @test dim(P) == 3
  @test length(vertices(P)) == 4

  res = s_characters(tbl, rational = false, ppow_nonzero = false);
  @test length(res) == 24
  res = s_characters(tbl, rational = true, ppow_nonzero = false);
  @test length(res) == 16
  res = s_characters(tbl, rational = false, ppow_nonzero = true);
  @test length(res) == 0
  res = s_characters(tbl, rational = true, ppow_nonzero = true);
  @test length(res) == 0

  tbl = character_table("2.A5");
  res = s_characters_via_factor_group(tbl, [1, 2], rational = false);
  @test length(res) == 106
  res = s_characters_via_factor_group(tbl, [1, 2], rational = true, ppow_nonzero = false);
  @test length(res) == 74
  res = s_characters_via_factor_group(tbl, [1, 2], rational = false, ppow_nonzero = true);
  @test length(res) == 0
  res = s_characters_via_factor_group(tbl, [1, 2], rational = true, ppow_nonzero = true);
  @test length(res) == 0
end
