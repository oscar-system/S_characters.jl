# utility function:
# save the result of `s_characters` to a file
function save_s_characters(filename::String, s_chars::Vector)
  l = [coordinates(ZZRingElem, v) for v in s_chars]
  tbl = parent(s_chars[1])
  info = Dict(:name => identifier(tbl), :coordinates => l)
  save(filename, info)
end


# utility function:
# create the LaTeX table in the paper from the stored data
names = [
          [ # alternating groups and their extensions
            [ "A8", "\\fA_{8}" ],
            [ "2.A8", "2.\\fA_{8}" ],
            [ "A9", "\\fA_{9}" ],
            [ "2.A9", "2.\\fA_{9}" ],
            [ "A10", "\\fA_{10}" ],
            [ "2.A10", "2.\\fA_{10}" ],
            [ "A11", "\\fA_{11}" ],
            [ "2.A11", "2.\\fA_{11}" ],
          ],
          [ # sporadic simple groups and their extensions
            [ "M12", "M_{12}" ],
            [ "2.M12", "2.M_{12}" ],
            [ "M12.2", "M_{12}.2" ],
            [ "J1", "J_{1}" ],
            [ "J2", "J_{2}" ],
            [ "2.J2", "2.J_{2}" ],
            [ "J2.2", "J_{2}.2" ],
            [ "J3", "J_{3}" ],
            [ "McL", "M^{c}L" ],
            [ "HS", "HS" ],
            [ "2.HS", "2.HS" ],
            [ "M24", "M_{24}" ],
          ],
          [ # simple groups of Lie type and their extensions
            [ "L4(3)", "\\PSL_{4}(3)" ],
            [ "L5(2)", "\\PSL_{5}(2)" ],
            [ "L5(2).2", "\\PSL_{5}(2).2" ],
            [ "S4(4)", "\\PSp_{4}(4)" ],
            [ "S4(4).2", "\\PSp_{4}(4).2" ],
            [ "S4(5)", "\\PSp_{4}(5)" ],
            [ "S6(2)", "\\PSp_{6}(2)" ],
            [ "U4(2)", "\\PSU_{4}(2)" ],
            [ "2.U4(2)", "2.\\PSU_{4}(2)" ],
            [ "2F4(2)'", "{}^{2}F_{4}(2)'" ],
            [ "2F4(2)'.2", "{}^{2}F_{4}(2)'.2" ],
            [ "R(27)", "^2G_2(27)" ],
            [ "G2(3)", "G_{2}(3)" ],
          ],
          [ # some other perfect groups
            [ "2^4:a8", "2^{4}\\!:\\!\\fA_{8}" ],
            [ "2^4.a8", "2^{4}.\\fA_{8}" ],
            [ "2^5.psl(5,2)", "2^{5}.\\PSL_{5}(2)" ],
            [ "2^6:A8", "2^{6}\\!:\\!\\fA_{8}" ],
          ],
        ];

function s_characters_table(names)
  # table header
  println("\\begin{table}")
  println("\\caption{Examples}   \\label{tab:ex}")
  println("\\(\n\\begin{array}{lrrrrr}\n  \\toprule")
  println("  G & \\# \\text{classes} & \\# \\text{real} & \\# \\text{rat.} & ",
          "\\# \\text{S-char.}\n",
          "  & \\# \\text{virt.\\ S-char.} \\\\\n",
          "  \\midrule")

  # run over the name portions
  for l in names
    for pair in l
      name = pair[1]
      latexname = pair[2]

      r = load(joinpath(@__DIR__, "../data", name))
      list = r[:coordinates]
      tbl = character_table(name)

      # compute the *faithful* example characters,
      # i.e. those which do not live in a proper factor group.
      # Note that we cannot use `class_positions_of_kernel` for
      # virtual characters.
      nsg = setdiff(class_positions_of_normal_subgroups(tbl), [[1]])
      nccl = length(tbl)
      faithirrpos = [filter(i -> !is_subset(n, class_positions_of_kernel(tbl[i])), 1:nccl) for n in nsg]
      faith = filter(x -> all(l -> any(i -> x[i] != 0, l), faithirrpos), list)

      # compute the *faithful* ordinary ones:
      ord = filter(x -> minimum(x) >= 0, list)
      ordfaith = intersect(ord, faith)

      # compute the values for the first three columns
      nreal = length(unique([x + conj(x) for x in tbl]))
      nrat = length(unique(map(galois_orbit_sum, collect(tbl))))

      # print the table row for this group,
      # separate portions
      if length(faith) == 0
        print("% ")
      else
        print("  ")
      end
      println(latexname, " & ", nccl, " & ", nreal, " & ", nrat, " & ",
              length(faith), " & ", length(faith) - length(ordfaith), " \\\\",
              (pair == l[end] && l != names[end]) ? "%\n  [1.6ex]\n  %" : "" )
    end
  end

  # table footer
  println("  \\bottomrule\n\\end{array}\n\\)\n\\end{table}")
end
