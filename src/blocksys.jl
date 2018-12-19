#= Copyright: Konrad Kleczkowski =#

module blocksys
  using SparseArrays, LinearAlgebra
  
  ### Wczytuje, zgodnie ze specyfikacją, plik z macierzą.
  ###
  ### Zwraca trójkę, która zawiera kolejno macierz,
  ### wymiar macierzy i wymiar podmacierzy
  function load_matrix(filename :: String) 
    f = open(filename, "r")

    first_line = readline(f)

    n, l = split(first_line, " ")
    n = parse(Int64, n)
    l = parse(Int64, l)

    matA = spzeros(Float64, n, n)

    for line in eachline(f) 
      i, j, value = split(line, " ")
      i = parse(Int64, i)
      j = parse(Int64, j)
      value = parse(Float64, value)
      matA[i, j] = value
    end

    close(f)
    return (matA, n, l)
  end

  ### Wczytuje, zgodnie ze specyfikacją, 
  ### plik z wektorem prawych stron.
  ###
  ### Zwraca trójkę, która zawiera wektor 
  ### i wymiar wektora.
  function load_b(filename :: String)
    f = open(filename, "r")
    
    n = parse(Int64, readline(f))
    b = zeros(n)
    for i = 1 : n 
      line = readline(f)
      b[i] = parse(Float64, line)
    end
    close(f)
    return (b, n)
  end

  ### Wyisuje treść wektora do pliku tekstowego.
  ###
  ### Przyjmuje wektor do wypisania, nazwę pliku
  ### i zmienną boolowską opisującą, czy wypisać
  ### błąd względem wektora jednostkowego.
  function write_x(x :: Vector{Float64}, filename :: String, rel_err :: Bool = false)
    f = open(filename, "w")
    if rel_err
      o = ones(length(x))
      e = norm(x - o) / norm(o)
      println(f, e)
    end
    for xi in x 
      println(f, xi)
    end
    close(f)
  end

  ### Generuje wektor b na podstawie wczytanej macierzy A.
  function generate_b(A, n :: Int, l :: Int)
    b = zeros(n)
    
    # pierwszy blok
    for i = 1 : l
      s = 0.0
      for j = 1 : l
        s += A[i, j]
      end
      s += A[i, i + l]
      b[i] =s
    end

    # środkowe bloki
    for k = 2 : Int64(n / l) - 1
      for i = (k - 1) * l + 1 : k * l
        s = 0.0
        for j = (k - 1) * l : k * l
          s += A[i, j]
        end
        s += A[i, i + l]
        b[i] = s
      end
    end

    # ostatni blok
    for i = n - l + 1 : n
      s = 0.0
      for j = n - l : n
        s += A[i, j]
      end
      b[i] = s
    end

    return b
  end

  ### Wykonuje eliminację Gaussa, zwracając 
  ### macierz dolno trójkątną i modyfikując A
  ### tak, by otrzymać macierz górno trójkątną.
  function gauss(A, b :: Vector{Float64}, n :: Int, l :: Int)
    L = spzeros(Float64, n, n)
    for k = 1 : n - 1
      L[k, k] = 1.0
      for i = k + 1 : k + l - k % l
        L[i, k] = A[i, k] / A[k, k]
        A[i, k] = 0.0
        for j = k + 1 : min(n, k + 2 * l)
          A[i, j] = A[i, j] - L[i, k] * A[k, j]
        end
        b[i] = b[i] - L[i, k] * b[k]
      end
    end
    L[n, n] = 1.0
    return L
  end
  
  ### Rozwiązuje układ z macierzą górno trójkątną.
  function solve_ut_matrix(U, b :: Vector{Float64}, n :: Int, l :: Int)
    x = zeros(n)
    for i = n : -1 : 1
      s = 0.0
      for j = i + 1 : min(n, i + 2 * l)
        s = s + U[i, j] * x[j]
      end
      x[i] = (b[i] - s) / U[i, i]
    end
    return x
  end

  ### Rozwiązuje układ z macierzą dolno trójkątną.
  function solve_lt_matrix(L, b :: Vector{Float64}, n :: Int, l :: Int)
    x = zeros(n)
    count = 0
    shift = 0
    for i = 1 : n
      s = 0.0 
      
      for j = max(1, trunc(Int64, (i - 1) / l) * l) : i - 1
        s = s + L[i, j] * x[j]
      end

      x[i] = (b[i] - s) / L[i, i]

      if count % l == 0 
        shift += l
      end
      count += 1
    end
    return x
  end
  
  ### Wykonuje eliminację Gaussa z częściowym wyborem.
  function choose_gauss(A, b :: Vector{Float64}, 
                        n :: Int, l :: Int, eps :: Float64)
    sigma = zeros(Int64, n) 
    for k = 1 : n
      sigma[k] = k
    end
    L = spzeros(Float64, n, n)
    for k = 1 : n - 1
      max = abs(A[k, k])
      max_row = k
      for i = k + 1 : k + l - k % l
        if abs(A[i, k]) > max 
          max = abs(A[i, k])
          max_row = i
        end
      end
      if max_row > k
        temp = spzeros(Float64, n)
        for t = k : min(n, k + 2 * l)
          temp[t] = A[k, t]
        end
        for t = k : min(n, k + 2 * l)
          A[k, t] = A[max_row, t]
        end
        for t = k : min(n, k + 2 * l)
          A[max_row, t] = temp[t]
        end
        for t =  Base.max(1, trunc(Int64, (k - 1) / l) * l) : k 
          temp[t] = L[k, t]
        end
        for t = Base.max(1, trunc(Int64, (k - 1) / l) * l) : k 
          L[k, t] = L[max_row, t]
        end
        for t = Base.max(1, trunc(Int64, (k - 1) / l) * l) : k
          L[max_row, t] = temp[t] 
        end

        b[max_row], b[k] = b[k], b[max_row]
        sigma[max_row], sigma[k] = sigma[k], sigma[max_row]
      end
      L[k, k] = 1.0
      for i = k + 1: k + l - k % l
        L[i, k] = A[i, k] / A[k, k]
        A[i, k] = 0.0
        for j = k + 1 : min(n, k + 2 * l)
          A[i, j] = A[i, j] - L[i, k] * A[k, j]
        end
        b[i] = b[i] - L[i, k] * b[k]
      end
    end
    L[n, n] = 1.0
    return (L, sigma)
  end

  ### Zwraca macierz A rozłożoną na macierze L i U.
  function compute_lu(A,
                      n :: Int, l :: Int)
    U = copy(A)
    L = gauss(U, ones(n), n, l)
    return (L, U)
  end

  ### Zwraca macierz A rozłożoną na macierze L i U
  ### za pomocą eliminacji Gaussa z częściowym wyborem.
  function compute_choose_lu(A, 
                             n :: Int, l :: Int, eps :: Float64)
    U = copy(A)
    b = ones(n)
    L, sigma = choose_gauss(U, b, n, l, eps)
    return (L, U, sigma)
  end

  ### Rozwiązuje układ równań bezpośrednio 
  ### z algorytmu eliminacji Gaussa.
  function solve_gauss(A, b :: Vector{Float64},
                       n :: Int, l :: Int)
    U = copy(A)
    bb = copy(b)
    L = gauss(U, bb, n, l)
    x = solve_ut_matrix(U, bb, n, l)
    return x
  end

  ### Rozwiązuje układ równań bezpośrednio
  ### z algorytmu eliminacji Gaussa
  ### z częściowym wyborem.
  function solve_choose_gauss(A, b :: Vector{Float64},
                              n :: Int, l :: Int, eps :: Float64)
    U = copy(A)
    bp = copy(b)
    L, sigma = choose_gauss(U, bp, n, l, eps)
    x = solve_ut_matrix(U, bp, n, l)

    return x
  end

  ### Rozwiązuje układ równań za pomocą rozłożonej macierzy
  ### do L i U.
  function solve_lu(L, U, b :: Vector{Float64}, 
                    n :: Int, l :: Int)
    y = solve_lt_matrix(L, b, n, l)
    x = solve_ut_matrix(U, y, n, l)
    return x
  end

  ### Rozwiązuje układ równań za pomocą rozłożonej macierzy
  ### do L i U i odwraca spermutowany wynikowy wektor
  ### na podstawie wygenerowanej permutacji sigma.
  function solve_choose_lu(L, U, sigma, b :: Vector{Float64}, 
                           n :: Int, l :: Int)
    bp = zeros(n)
    for i = 1 : n 
      bp[i] = b[sigma[i]]
    end
    y = solve_lt_matrix(L, bp, n, l)
    x = solve_ut_matrix(U, y, n, l)
    return x
  end
end
