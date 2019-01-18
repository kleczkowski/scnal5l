#= Copyright: Konrad Kleczkowski =#

include("blocksys.jl")

using LinearAlgebra
using SparseArrays

### Uruchamia test dla:
###   * eliminacji Gaussa bez i z częściowym wyborem (eps = 0.1)
###   * rozwiązywania układu równań na podstawie rozkładu LU
function run_basic_solve_test(mat_file :: String,
                              b_file :: String,
                              out_file :: String)
  println("--- $(out_file):")
  A, n, l = blocksys.load_matrix(mat_file)
  b, _ = blocksys.load_b(b_file)
  println("solve_gauss:")
  @time xg = blocksys.solve_gauss(A, b, n, l)
  println("solve_choose_gauss:")
  @time xcg = blocksys.solve_choose_gauss(A, b, n, l, 0.1)

  blocksys.write_x(xg, "$(out_file)_gauss.txt", true)
  blocksys.write_x(xcg, "$(out_file)_choose_gauss.txt", true)

  println("compute_lu:")
  @time Lb, Ub = blocksys.compute_lu(A, n, l)
  println("  błąd względny: $(norm(A - Lb * Ub) / norm(A))")
  println("solve_lu:")
  @time xblu = blocksys.solve_lu(Lb, Ub, b, n, l)

  println("compute_choose_lu:")
  @time Lc, Uc, sigma = blocksys.compute_choose_lu(A, n, l, 0.1)
  id = zeros(Int64, n)
  for i = 1 : n
    id[i] = i
  end
  Ap = permute(A, sigma, id)
  println("  błąd względny: $(norm(Ap - Lc * Uc) / norm(Ap))")
  println("solve_choose_lu:")
  @time xclu = blocksys.solve_choose_lu(Lc, Uc, sigma, b, n, l)
  blocksys.write_x(xg, "$(out_file)_lu.txt", true)
  blocksys.write_x(xcg, "$(out_file)_choose_lu.txt", true)
end

#run_basic_solve_test("../resources/Dane16_1_1/A.txt",
#                     "../resources/Dane16_1_1/b.txt",
#                     "dane16")
#run_basic_solve_test("../resources/Dane10000_1_1/A.txt",
#                     "../resources/Dane10000_1_1/b.txt",
#                     "dane10000")
#run_basic_solve_test("../resources/Dane50000_1_1/A.txt",
#                     "../resources/Dane50000_1_1/b.txt",
#                     "dane50000")
