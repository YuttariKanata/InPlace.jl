using InPlace
using BenchmarkTools
#using LinearAlgebra
using Base.Threads

# --- [1. Allocation & Throughput: InPlace vs Standard] ---

function benchmark_basic_arithmetic(N = 1_000_000)
    println("--- Test 1: Micro-benchmark (Arithmetic) ---")
    GC.gc()
    
    # 準備
    a, b = BigInt(rand(UInt128)), BigInt(rand(UInt128))
    res = BigInt(0)

    println("Standard Julia (res = a + b):")
    GC.gc()
    
    # 標準演算は毎回新しいオブジェクトを生成するためアロケーションが発生
    @btime for i in 1:$N; $res = $a + $b; end

    println("InPlace.jl (add!(res, a, b)):")
    # InPlaceは既存のバッファを再利用するためアロケーションは0
    @btime for i in 1:$N; add!($res, $a, $b); end
    println()
end

# --- [2. Memory Pressure & GC Impact Test] ---
# メモリの少ない環境や内蔵GPU環境で重要な「GCの介入」を可視化

function benchmark_gc_pressure(N = 2_000_000)
    println("--- Test 2: Memory Pressure & GC Impact ---")
    GC.gc()
    
    # 標準演算での実行
    t_std = @elapsed begin
        s = BigInt(0)
        for i in 1:N
            s = s + i  # 毎回一時的なBigIntが生成される
        end
    end
    println("Standard Total Time: ", round(t_std, digits=3), "s")
    GC.gc()

    # InPlaceでの実行
    t_ip = @elapsed begin
        s_ip = BigInt(0)
        for i in 1:N
            add!(s_ip, s_ip, i) # 内部で _ui 系の最適化パスが通る
        end
    end
    println("InPlace Total Time:  ", round(t_ip, digits=3), "s")
    println("Speedup: ", round(t_std / t_ip, digits=2), "x")
    println()
end

# --- [3. Multi-threaded Stress (CPU Parallelism)] ---
# 並列実行時にGCがボトルネックにならないかを検証

function benchmark_multithreaded(N = 500_000)
    println("--- Test 3: Multi-threaded Scaling (Threads: $(nthreads())) ---")
    GC.gc()

    # 標準演算（スレッドごとにGCが走る可能性がある）
    t_std = @elapsed begin
        @threads for t in 1:nthreads()
            acc = BigInt(0)
            for i in 1:N
                acc = acc + i
            end
        end
    end

    GC.gc()

    # InPlace（アロケーションがないためスレッド間の干渉が極小）
    t_ip = @elapsed begin
        @threads for t in 1:nthreads()
            acc = BigInt(0)
            for i in 1:N
                add!(acc, acc, i)
            end
        end
    end
    
    println("Standard Multi-thread: ", round(t_std, digits=3), "s")
    println("InPlace  Multi-thread: ", round(t_ip, digits=3), "s")
    println()
end

# --- [4. Precision Scaling (BigFloat)] ---
# 精度（bit数）を上げた際のスループットを計測

function benchmark_float_precision()
    println("--- Test 4: BigFloat Precision Scaling ---")
    
    for prec in [128, 512, 2048]
        GC.gc()
        println("Precision: $prec bits")
        x = BigFloat(1.1, precision=prec)
        y = BigFloat(2.2, precision=prec)
        r = BigFloat(0, precision=prec)
        
        # mpfr_mul を直接叩くことによる効率化を検証[cite: 6]
        print("  InPlace mul!: ")
        @btime mul!($r, $x, $y)
    end
    println()
end

# --- [Main Execution] ---

function run_all_tests()
    println("InPlace.jl General Performance Test Suite")
    println("="^40)
    
    benchmark_basic_arithmetic()
    benchmark_gc_pressure()
    benchmark_multithreaded()
    benchmark_float_precision()
    
    println("="^40)
    println("Test completed.")
end

run_all_tests()