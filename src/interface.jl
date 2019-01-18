#= Copyright: Konrad Kleczkowski =#

include("blocksys.jl")
include("test.jl")

function load_matrix_and_generate_b(inFile :: String, outFile :: String)
  A, n, l = blocksys.load_matrix(inFile)
  b = blocksys.generate_b(A, n, l)
  fb = open(outFile, "w")
  write(fb, "$(n)\n")
  for bi in b
    write(fb, "$(bi)\n")
  end
  close(fb)
end

mode = ARGS[1]

if mode == "--generate-b"
  inFile = ARGS[2]
  outFile = ARGS[3]
  load_matrix_and_generate_b(inFile, outFile)
end

if mode == "--run-test"
  matrixFile = ARGS[2]
  bFile = ARGS[3]
  outFile = ARGS[4]
  run_basic_solve_test(matrixFile, bFile, outFile)
end
