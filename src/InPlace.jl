module InPlace


include("gmp_ffi.jl")    # 純粋な ccall ラッパー
include("mpfr_ffi.jl")   # MPFR版 ccall ラッパー
include("arithmetic.jl") # add! などの破壊的インターフェース


end
