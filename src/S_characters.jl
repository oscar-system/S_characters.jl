module S_characters

using Oscar

export s_character_simplex
export s_characters
export s_characters_via_factor_group
export one_rational_s_character_via_milp

# utility function:
# For a vector `A` of vectors (class functions),
# compute a vector for which duplicate columns are removed,
# and return `M, fus, proj`
# where `M` is the vector of collapsed vectors,
# and `fus` and `proj` are vectors of integers
# that describe the relation between `A` and `M`.
function _collapse_columns(A::Vector)
  m = length(A)
  n = length(A[1])
  fus = collect(1:n)
  proj = Int[1]
  take = ones(Bool, n)
  for j in 2:n
    # find the first column equal to the j-th column
    for k in 1:(j-1)
      if all(v -> v[k] == v[j], A)
        fus[j] = fus[k]
        take[j] = false
        continue
      end
    end
    if take[j]
      push!(proj, j)
      fus[j] = length(proj)
    end
  end

  return map(v -> v[findall(take)], A), fus, proj
end

"""
    s_character_simplex(tbl::Oscar.GAPGroupCharacterTable;
               rational::Bool = true,
               ppow_nonzero::Bool = false,
               irrats::Symbol = :nf)

Return `P, galoissums, ppow_pos` where

- `P` is the polyhedron that is defined by the inequalities given by
  the rational irreducible characters of `tbl` (if `rational` is `true`)
  or the real irreducible characters of `tbl` (if `rational` is `false`),

- `galoissums` is the vector of real or rational irreducibles of `tbl`, and

- `ppow_pos` is a vector of column positions in `tbl` for which a strictly
  positive value is prescribed in the S-characters of `tbl` encoded by `P`;
  these are the rational classes of elements of prime power order if
  `ppow_nonzero` is `true`, otherwise `ppow_pos` is empty.

The optional argument `irrats` defines how the irrational entries of the
defining matrix of `P` are represented:
The value `:nf` means that these entries belong to a common embedded
number field (constructed from the common field of character values),
other values mean that these entries lie in `algebraic_closure(QQ)`.
If the degree of the common number field is small then the computations
are expected to be faster when `:nf` is chosen.

The polytope `P` is known to be a simplex.

# Examples
```jldoctest
julia> tbl = character_table("A5");

julia> P, galoissums, ppow_pos = s_character_simplex(tbl);

julia> P
Polyhedron in ambient dimension 3

julia> vertices(P)
4-element SubObjectIterator{PointVector{QQFieldElem}}:
 [3, 4, 5]
 [1//2, -1, 0]
 [-1, 0, 1]
 [0, 1, -1]
```
"""
function s_character_simplex(tbl::Oscar.GAPGroupCharacterTable;
                    rational::Bool = true,
                    ppow_nonzero::Bool = false,
                    irrats::Symbol = :nf)
  @req order(tbl) > 1 "the underlying group must not be trivial"

  # Compute the orbit sums on the nontrivial irreducibles
  # under complex conjugation if `rational == false`
  # and under all Galois automorphisms if `rational == true`;
  # omit the trivial character.
  triv = trivial_character(tbl)
  pos = findfirst(isequal(triv), tbl)
  if rational
    c = collect(tbl)
    nontriv = vcat(c[1:pos-1], c[pos+1:end])
    galoissums = unique(map(galois_orbit_sum, nontriv))
  else
    galoissums = Oscar.GAPGroupClassFunction[]
    for i in 1:length(tbl)
      i == pos && continue
      chi = tbl[i]
      psi = conj(chi)
      if psi == chi
        push!(galoissums, chi)
      elseif findfirst(isequal(psi), tbl) > i
        push!(galoissums, chi + psi)
      end
    end
  end

  # Collapse equal columns (do not duplicate the defining inequalities).
  coll, fus, proj = _collapse_columns(map(values, galoissums))

  # Negate and transpose the matrix,
  # map irrational entries to the field in question.
  m = length(proj)
  n = length(coll)
  ratcol = ones(Bool, m)
  if rational
    F = QQ
    mp = x -> F(x)
  elseif irrats == :nf
    CF, emb = character_field(galoissums)
    e = real_embeddings(CF)[1]
    F = Hecke.embedded_field(CF, e)[1]
    mp = x -> F(preimage(emb, x))
  else
    F = algebraic_closure(QQ)
    mp = x -> F(x)
  end
  T = elem_type(typeof(F))
  A = Matrix{T}(undef, m, n)
  for i in 1:m, j in 1:n
    x = coll[j][i]
    if x.c != 1
      ratcol[i] = false
    end
    A[i,j] = - mp(x)
  end

  # Initialize the right hand side `b` of the system of inequalities.
  b = ones(Int, m)
  if ppow_nonzero
    # Replace the `1` in `b` by `0` in all those positions of columns
    # that are rational and belong to elements of prime power order.
    ppow_pos = findall(x -> x == 1 || is_prime_power_with_data(x)[1],
                       orders_class_representatives(tbl))
    for i in unique(fus[ppow_pos])
      if rational || ratcol[i]
        b[i] = 0
      end
    end
  else
    # Do not prescribe a stronger condition for elements of prime power order.
    ppow_pos = Int[]
  end

  # Compute the polyhedron that is defined by the inequalities
  # given by the (real or rational) irreducible characters.
  # Its lattice points are the coefficient vectors of
  # all integer linear combinations `psi` of the real or rational irreducibles
  # such that `psi + 1` has nonnegative values
  # (and if applicable then `psi + 1` is nonzero on certain classes,
  # see above).
  P = polyhedron(F, A, b)

  return P, galoissums, ppow_pos
end


"""
    s_characters(tbl::Oscar.GAPGroupCharacterTable;
                 rational::Bool = true,
                 ppow_nonzero::Bool = false,
                 irrats::Symbol = :nf)

Return the vector of nontrivial S-characters of `tbl`
that are described by the arguments.
See [`s_character_simplex`](@ref) for the meaning of the keyword arguments.

# Examples
```jldoctest
julia> tbl = character_table("A5");

julia> res = s_characters(tbl, rational = false, ppow_nonzero = false);

julia> length(res)  # all nontrivial S-characters of A5
24

julia> res = s_characters(tbl, rational = true, ppow_nonzero = false);

julia> length(res)  # all rational nontrivial S-characters of A5
16

julia> res = s_characters(tbl, rational = false, ppow_nonzero = true);

julia> length(res)  # no nontrivial S-characters positive on p-power elements
0

julia> res = s_characters(tbl, rational = true, ppow_nonzero = true);

julia> length(res)  # of course no rational such S-characters
0
```
"""
function s_characters(tbl::Oscar.GAPGroupCharacterTable;
                       rational::Bool = true,
                       ppow_nonzero::Bool = false,
                       irrats::Symbol = :nf)
  # The trivial group has no nontrivial S-characters.
  order(tbl) == 1 && return typeof(tbl[1])[]

  P, galoissums, ppow_pos = s_character_simplex(tbl,
                              rational = rational,
                              ppow_nonzero = ppow_nonzero,
                              irrats = irrats)
  ll = lattice_points(P)

  # Compute the corresponding S-characters;
  # if `ppow_nonzero` is true then collect those that are nonzero
  # on all classes of elements of prime power order.
  res = typeof(galoissums[1])[]
  triv = trivial_character(parent(galoissums[1]))
  n = length(galoissums)
  for v in ll
    is_zero(v) && continue
    cand = triv + sum(v[i] * galoissums[i] for i in 1:n)
    if !ppow_nonzero || all(x -> !is_zero(cand[x]), ppow_pos)
      push!(res, cand)
    end
  end

  return res
end

@doc raw"""
    s_characters_via_factor_group(tbl::Oscar.GAPGroupCharacterTable,
        nclasses::Vector{Int};
        rational::Bool = true,
        ppow_nonzero::Bool = false,
        irrats::Symbol = :nf)

Return the vector of nontrivial S-characters of `tbl`
that are described by the arguments.

`nclasses` must a list of class positions in `tbl` that form a normal
subgroup $N$.
The function proceeds in two steps.
First it computes the S-characters in question for the factor group modulo $N$,
then it computes, for each such S-character $\psi$, those S-characters of `tbl`
whose projection to the factor group is $\psi$.

See [`s_character_simplex`](@ref) for the meaning of the keyword arguments.

# Examples
```jldoctest
julia> tbl = character_table("2.A5");

julia> res = s_characters_via_factor_group(tbl, [1, 2], rational = false);

julia> length(res)  # all nontrivial S-characters of 2.A5
106

julia> res = s_characters_via_factor_group(tbl, [1, 2], rational = true, ppow_nonzero = false);

julia> length(res)  # all rational nontrivial S-characters of A5
74

julia> res = s_characters_via_factor_group(tbl, [1, 2], rational = false, ppow_nonzero = true);

julia> length(res)  # no nontrivial S-characters positive on p-power elements
0

julia> res = s_characters_via_factor_group(tbl, [1, 2], rational = true, ppow_nonzero = true);

julia> length(res)  # of course no rational such S-characters
0
```
"""
function s_characters_via_factor_group(tbl::Oscar.GAPGroupCharacterTable,
             nclasses::Vector{Int};
             rational::Bool = true,
             ppow_nonzero::Bool = false,
             irrats::Symbol = :nf)

  # Compute the S-characters for the factor group in question.
  facttbl, factfus = quo(tbl, nclasses)
  fact_s_chars = s_characters(facttbl,
                       rational = rational,
                       ppow_nonzero = ppow_nonzero,
                       irrats = irrats)

  # Compute the orbit sums on the N-faithful irreducibles
  # under complex conjugation if `rational == false`
  # and under all Galois automorphisms if `rational == true`.
  if rational
    N_faith = filter(x -> !is_subset(nclasses, class_positions_of_kernel(x)),
                     collect(tbl))
    galoissums = unique(map(galois_orbit_sum, N_faith))
  else
    galoissums = Oscar.GAPGroupClassFunction[]
    for i in 1:length(tbl)
      chi = tbl[i]
      is_subset(nclasses, class_positions_of_kernel(chi)) && continue
      psi = conj(chi)
      if psi == chi
        push!(galoissums, chi)
      else
        pos = findfirst(isequal(psi), tbl)
        if pos > i
          push!(galoissums, chi + psi)
        end
      end
    end
  end

  # Collapse equal columns (do not duplicate the defining inequalities).
  coll, fus, proj = _collapse_columns(map(values, galoissums))

  # Negate and transpose the matrix,
  # map irrational entries to the field in question.
  m = length(proj)
  n = length(coll)
  ratcol = ones(Bool, m)
  if rational
    F = QQ
    mp = x -> F(x)
  elseif irrats == :nf
    # We may need a proper field extension.
    # Make sure that the field contains also
    # the character values for the factor group.
    # (Since we consider only *one* solution for the factor group
    # at a time, it is only necessary to consider the field of values
    # of this vector.
    # All solutions we found up to now are rational,
    # thus we need not think about such an improved strategy yet.)
    CF, emb = character_field(vcat(galoissums, fact_s_chars))
    e = real_embeddings(CF)[1]
    F = Hecke.embedded_field(CF, e)[1]
    mp = x -> F(preimage(emb, x))
  else
    F = algebraic_closure(QQ)
    mp = x -> F(x)
  end
  T = elem_type(typeof(F))
  A = Matrix{T}(undef, m, n)
  for i in 1:m, j in 1:n
    x = coll[j][i]
    if x.c != 1
      ratcol[i] = false
    end
    A[i,j] = - mp(x)
  end

  if ppow_nonzero
    # Decrease the right hand side by `1` in all those positions of columns
    # that are rational and belong to elements of prime power order.
    ppow_pos = findall(x -> x == 1 || is_prime_power_with_data(x)[1],
                       orders_class_representatives(tbl))
  else
    # Do not prescribe a stronger condition for elements of prime power order.
    ppow_pos = Int[]
  end

  # For each solution v, compute the solutions in the next smaller
  # factor group that lie above v.
  res = [restrict(chi, tbl) for chi in fact_s_chars]

  # Solutions for `tbl` may project to the trivial character of the factor.
  pushfirst!(fact_s_chars, trivial_character(facttbl))

  for v in fact_s_chars
    # Initialize the right hand side `b` of the system of inequalities.
    # (Fill up v by characters not from the factor group)
    # First inflate v, then collapse it.
# Improve this:
# Always take the strongest condition from the set of columns that map to the
# same column of the factor group.
    infl_v = values(v)[factfus]
    coll_v = infl_v[proj]
    b = map(mp, coll_v)
    if ppow_nonzero
      for i in unique(fus[ppow_pos])
        if rational || ratcol[i]
          b[i] = b[i] - 1
        end
      end
    end

    P = polyhedron(F, A, b)
    ll = lattice_points(P)

    # Compute the corresponding S-characters;
    # if `ppow_nonzero` is `true` then collect those that are nonzero
    # on all classes of elements of prime power order.
    factcand = Oscar.class_function(tbl, infl_v)
    n = length(galoissums)
    for w in ll
      is_zero(w) && continue
      cand = sum(w[i] * galoissums[i] for i in 1:n; init = factcand)
      if !ppow_nonzero || all(x -> !is_zero(cand[x]), ppow_pos)
        push!(res, cand)
      end
    end
  end

  return res
end;

@doc raw"""
    one_rational_s_character_via_milp(tbl::Oscar.GAPGroupCharacterTable)

Return `flag, psi` where `flag` is `true` if the program has found
a nontrivial S-character of `tbl` that is nonzero on all elements of
prime power order, and `false` otherwise.

If `flag` is `true` then `psi` is the S-character that was found,
otherwise `psi` is the trivial character of `tbl`.

# Examples
```jldoctest
julia> flag, psi = one_rational_s_character_via_milp(character_table("A5"));

julia> flag
false

julia> flag, psi = one_rational_s_character_via_milp(character_table("A8"));

julia> flag
true

julia> println(coordinates(psi))
QQFieldElem[1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3]
```
"""
function one_rational_s_character_via_milp(tbl::Oscar.GAPGroupCharacterTable)
  P, galoissums, ppow_pos = s_character_simplex(tbl,
                              rational = true,
                              ppow_nonzero = true,
                              irrats = :nf);
  d = length(galoissums)
  milp = mixed_integer_linear_program(P, repeat([1], d))
  _, sol = solve_milp(milp)
  if is_zero(sol)
    return false, trivial_character(tbl)
  else
    # By construction, we know that the S-character is positive on `ppow_pos`.
    cand = sum(ZZ(sol[i]) * galoissums[i] for i in 1:d; init = trivial_character(tbl))
    return true, cand
  end
end

end # module S_characters
