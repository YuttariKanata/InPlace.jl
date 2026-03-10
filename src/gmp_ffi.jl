# src/gmp_ffi.jl
# Low-level GMP bindings using GMP_jll

using GMP_jll
using Base.GMP.MPZ: mpz_t, bitcnt_t
using Base.GMP: CulongMax, ClongMax, CdoubleMax, Limb
using Libdl

const Cbitcnt = bitcnt_t
const mp_limb_t = Limb  # mp_limb_t の実体を確認 (通常は UInt64)
const mp_size_t = BigInt.types[2]



if Cint == Int32
    const CintMax = Union{Int8, Int16, Int32}
else
    const CintMax = Union{Int8, Int16, Int32, Int64}
end
if bitcnt_t == UInt32
    const CbitcntMax = Union{UInt8, UInt16, UInt32}
else
    const CbitcntMax = Union{UInt8, UInt16, UInt32, UInt64}
end
if Csize_t == UInt32
    const Csize_tMax = Union{UInt8, UInt16, UInt32}
else
    const Csize_tMax = Union{UInt8, UInt16, UInt32, UInt64}
end
if mp_size_t == Int32
    const Cmp_size_tMax = Union{Int8, Int16, Int32}
else
    const Cmp_size_tMax = Union{Int8, Int16, Int32, Int64}
end
if mp_limb_t == UInt32
    const Cmp_limb_tMax = Union{UInt8, UInt16, UInt32}
else
    const Cmp_limb_tMax = Union{UInt8, UInt16, UInt32, UInt64}
end



# --- Type Mapping for Arguments ---
const TYPE_MAP = Dict(
    :M   => (BigInt,        :mpz_t         ), # mpz_t
    :UI  => (CulongMax,     :Culong        ), # unsigned long
    :SI  => (ClongMax,      :Clong         ), # signed long
    :I   => (CintMax,       :Cint          ), # int
    :D   => (CdoubleMax,    :Cdouble       ), # double
    :BC  => (CbitcntMax,    :bitcnt_t      ), # bit count (Culong)
    :S   => (Cstring,       :Cstring       ), # char*
    :PSI => (Ref{Clong},    :(Ptr{Clong})  ), # signed long int *
    :PUI => (Ref{Culong},   :(Ptr{Culong}) ), # unsigned long int *
    :SZ  => (Csize_tMax,    :Csize_t       ), # size_t
    :ML  => (Cmp_limb_tMax, :mp_limb_t     ), # mp_limb_t (ほとんどの64bit環境でCsize_tと同等)
    :MS  => (Cmp_size_tMax, :mp_size_t     ), # mp_size_t
)



# --- GMP Function Definitions ---
# (戻り値の型, 関数名, 引数パターンのリスト)
const GMP_DEFINITIONS = [
    # 5.2 Assignment Functions                                     
    (:Cvoid,     :mpz_set,                 [:M,   :M              ]),
    (:Cvoid,     :mpz_set_ui,              [:M,   :UI             ]),
    (:Cvoid,     :mpz_set_si,              [:M,   :SI             ]),
    (:Cvoid,     :mpz_set_d,               [:M,   :D              ]),
    (:Cvoid,     :mpz_swap,                [:M,   :M              ]),
                                                                   
    # 5.4 Conversion Functions                                     
    (:Culong,    :mpz_get_ui,              [:M                    ]),
    (:Clong,     :mpz_get_si,              [:M                    ]),
    (:Cdouble,   :mpz_get_d,               [:M                    ]),
    (:Cdouble,   :mpz_get_d_2exp,          [:PSI, :M              ]),
                                                                   
    # 5.5 Arithmetic Functions                                     
    (:Cvoid,     :mpz_add,                 [:M,   :M,   :M        ]),
    (:Cvoid,     :mpz_add_ui,              [:M,   :M,   :UI       ]),
    (:Cvoid,     :mpz_sub,                 [:M,   :M,   :M        ]),
    (:Cvoid,     :mpz_sub_ui,              [:M,   :M,   :UI       ]),
    (:Cvoid,     :mpz_ui_sub,              [:M,   :UI,  :M        ]),
    (:Cvoid,     :mpz_mul,                 [:M,   :M,   :M        ]),
    (:Cvoid,     :mpz_mul_si,              [:M,   :M,   :SI       ]),
    (:Cvoid,     :mpz_mul_ui,              [:M,   :M,   :UI       ]),
    (:Cvoid,     :mpz_addmul,              [:M,   :M,   :M        ]),
    (:Cvoid,     :mpz_addmul_ui,           [:M,   :M,   :UI       ]),
    (:Cvoid,     :mpz_submul,              [:M,   :M,   :M        ]),
    (:Cvoid,     :mpz_submul_ui,           [:M,   :M,   :UI       ]),
    (:Cvoid,     :mpz_mul_2exp,            [:M,   :M,   :BC       ]),
    (:Cvoid,     :mpz_neg,                 [:M,   :M              ]),
    (:Cvoid,     :mpz_abs,                 [:M,   :M              ]),
                                                                   
    # 5.6 Division Functions                                       
    (:Cvoid,     :mpz_cdiv_q,              [:M,   :M,   :M        ]),
    (:Cvoid,     :mpz_cdiv_r,              [:M,   :M,   :M        ]),
    (:Cvoid,     :mpz_cdiv_qr,             [:M,   :M,   :M,   :M  ]),
    (:Culong,    :mpz_cdiv_q_ui,           [:M,   :M,   :UI       ]),
    (:Culong,    :mpz_cdiv_r_ui,           [:M,   :M,   :UI       ]),
    (:Culong,    :mpz_cdiv_qr_ui,          [:M,   :M,   :M,   :UI ]),
    (:Culong,    :mpz_cdiv_ui,             [:M,   :UI             ]),
    (:Cvoid,     :mpz_cdiv_q_2exp,         [:M,   :M,   :BC       ]),
    (:Cvoid,     :mpz_cdiv_r_2exp,         [:M,   :M,   :BC       ]),
    (:Cvoid,     :mpz_fdiv_q,              [:M,   :M,   :M        ]),
    (:Cvoid,     :mpz_fdiv_r,              [:M,   :M,   :M        ]),
    (:Cvoid,     :mpz_fdiv_qr,             [:M,   :M,   :M,   :M  ]),
    (:Culong,    :mpz_fdiv_q_ui,           [:M,   :M,   :UI       ]),
    (:Culong,    :mpz_fdiv_r_ui,           [:M,   :M,   :UI       ]),
    (:Culong,    :mpz_fdiv_qr_ui,          [:M,   :M,   :M,   :UI ]),
    (:Culong,    :mpz_fdiv_ui,             [:M,   :UI             ]),
    (:Cvoid,     :mpz_fdiv_q_2exp,         [:M,   :M,   :BC       ]),
    (:Cvoid,     :mpz_fdiv_r_2exp,         [:M,   :M,   :BC       ]),
    (:Cvoid,     :mpz_tdiv_q,              [:M,   :M,   :M        ]),
    (:Cvoid,     :mpz_tdiv_r,              [:M,   :M,   :M        ]),
    (:Cvoid,     :mpz_tdiv_qr,             [:M,   :M,   :M,   :M  ]),
    (:Culong,    :mpz_tdiv_q_ui,           [:M,   :M,   :UI       ]),
    (:Culong,    :mpz_tdiv_r_ui,           [:M,   :M,   :UI       ]),
    (:Culong,    :mpz_tdiv_qr_ui,          [:M,   :M,   :M,   :UI ]),
    (:Culong,    :mpz_tdiv_ui,             [:M,   :UI             ]),
    (:Cvoid,     :mpz_tdiv_q_2exp,         [:M,   :M,   :BC       ]),
    (:Cvoid,     :mpz_tdiv_r_2exp,         [:M,   :M,   :BC       ]),
    (:Cvoid,     :mpz_mod,                 [:M,   :M,   :M        ]),
    (:Culong,    :mpz_mod_ui,              [:M,   :M,   :UI       ]),
    (:Cvoid,     :mpz_divexact,            [:M,   :M,   :M        ]),
    (:Cvoid,     :mpz_divexact_ui,         [:M,   :M,   :UI       ]),
    (:Cint,      :mpz_divisible_p,         [:M,   :M              ]),
    (:Cint,      :mpz_divisible_ui_p,      [:M,   :UI             ]),
    (:Cint,      :mpz_divisible_2exp_p,    [:M,   :BC             ]),
    (:Cint,      :mpz_congruent_p,         [:M,   :M,   :M        ]),
    (:Cint,      :mpz_congruent_ui_p,      [:M,   :UI,  :UI       ]),
    (:Cint,      :mpz_congruent_2exp_p,    [:M,   :M,   :BC       ]),
                                                                   
    # 5.7 Exponentiation Functions                                 
    (:Cvoid,     :mpz_powm,                [:M,   :M,   :M,   :M  ]),
    (:Cvoid,     :mpz_powm_ui,             [:M,   :M,   :UI,  :M  ]),
    (:Cvoid,     :mpz_powm_sec,            [:M,   :M,   :M,   :M  ]),
    (:Cvoid,     :mpz_pow_ui,              [:M,   :M,   :UI       ]),
    (:Cvoid,     :mpz_ui_pow_ui,           [:M,   :UI,  :UI       ]),
                                                                   
    # 5.8 Root Extraction Functions                                
    (:Cint,      :mpz_root,                [:M,   :M,   :UI       ]),
    (:Cvoid,     :mpz_rootrem,             [:M,   :M,   :M,   :UI ]),
    (:Cvoid,     :mpz_sqrt,                [:M,   :M              ]),
    (:Cvoid,     :mpz_sqrtrem,             [:M,   :M,   :M        ]),
    (:Cint,      :mpz_perfect_power_p,     [:M,                   ]),
    (:Cint,      :mpz_perfect_square_p,    [:M,                   ]),
                                                                   
    # 5.9 Number Theoretic                                         
    (:Cint,      :mpz_probab_prime_p,      [:M,   :I              ]),
    (:Cvoid,     :mpz_nextprime,           [:M,   :M              ]),
    (:Cint,      :mpz_prevprime,           [:M,   :M              ]),
    (:Cvoid,     :mpz_gcd,                 [:M,   :M,   :M        ]),
    (:Culong,    :mpz_gcd_ui,              [:M,   :M,   :UI       ]),
    (:Cvoid,     :mpz_gcdext,              [:M,   :M,   :M,   :M  ]),
    (:Cvoid,     :mpz_lcm,                 [:M,   :M,   :M        ]),
    (:Cvoid,     :mpz_lcm_ui,              [:M,   :M,   :UI       ]),
    (:Cvoid,     :mpz_invert,              [:M,   :M,   :M        ]),
    (:Cint,      :mpz_jacobi,              [:M,   :M              ]),
    (:Cint,      :mpz_legendre,            [:M,   :M,   :M        ]),
    (:Cint,      :mpz_kronecker,           [:M,   :M              ]),
    (:Cint,      :mpz_kronecker_si,        [:M,   :SI             ]),
    (:Cint,      :mpz_kronecker_ui,        [:M,   :UI             ]),
    (:Cint,      :mpz_si_kronecker,        [:SI,  :M              ]),
    (:Cint,      :mpz_ui_kronecker,        [:UI,  :M              ]),
    (:Cbitcnt,   :mpz_remove,              [:M,   :M,   :M        ]),
    (:Cvoid,     :mpz_fac_ui,              [:M,   :UI,            ]),
    (:Cvoid,     :mpz_2fac_ui,             [:M,   :UI,            ]),
    (:Cvoid,     :mpz_mfac_uiui,           [:M,   :UI,  :UI       ]),
    (:Cvoid,     :mpz_primorial_ui,        [:M,   :UI             ]),
    (:Cvoid,     :mpz_bin_ui,              [:M,   :M,   :UI       ]),
    (:Cvoid,     :mpz_bin_uiui,            [:M,   :UI,  :UI       ]),
    (:Cvoid,     :mpz_fib_ui,              [:M,   :UI             ]),
    (:Cvoid,     :mpz_fib2_ui,             [:M,   :M,   :UI       ]),
    (:Cvoid,     :mpz_lucnum_ui,           [:M,   :UI             ]),
    (:Cvoid,     :mpz_lucnum2_ui,          [:M,   :M,   :UI       ]),
                                                                   
    # 5.10 Comparison Functions                                    
    (:Cint,      :mpz_cmp,                 [:M,   :M              ]),
    (:Cint,      :mpz_cmp_d,               [:M,   :D              ]),
    (:Cint,      :mpz_cmp_si,              [:M,   :SI             ]),
    (:Cint,      :mpz_cmp_ui,              [:M,   :UI             ]),
    (:Cint,      :mpz_cmpabs,              [:M,   :M              ]),
    (:Cint,      :mpz_cmpabs_d,            [:M,   :D              ]),
    (:Cint,      :mpz_cmpabs_ui,           [:M,   :UI             ]),
    (:Cint,      :mpz_sgn,                 [:M,                   ]),
                                                                   
    # 5.11 Logical and Bit Manipulation Functions                  
    (:Cvoid,     :mpz_and,                 [:M,   :M,   :M        ]),
    (:Cvoid,     :mpz_ior,                 [:M,   :M,   :M        ]),
    (:Cvoid,     :mpz_xor,                 [:M,   :M,   :M        ]),
    (:Cvoid,     :mpz_com,                 [:M,   :M              ]),
    (:Cbitcnt,   :mpz_popcount,            [:M,                   ]),
    (:Cbitcnt,   :mpz_hamdist,             [:M,   :M              ]),
    (:Cbitcnt,   :mpz_scan0,               [:M,   :BC             ]),
    (:Cbitcnt,   :mpz_scan1,               [:M,   :BC             ]),
    (:Cvoid,     :mpz_setbit,              [:M,   :BC             ]),
    (:Cvoid,     :mpz_clrbit,              [:M,   :BC             ]),
    (:Cvoid,     :mpz_combit,              [:M,   :BC             ]),
    (:Cint,      :mpz_tstbit,              [:M,   :BC             ]),
                                                                   
    # 5.15 Miscellaneous Function                                  
    (:Cint,      :mpz_fits_ulong_p,        [:M,                   ]),
    (:Cint,      :mpz_fits_slong_p,        [:M,                   ]),
    (:Cint,      :mpz_fits_uint_p,         [:M,                   ]),
    (:Cint,      :mpz_fits_sint_p,         [:M,                   ]),
    (:Cint,      :mpz_fits_ushort_p,       [:M,                   ]),
    (:Cint,      :mpz_fits_sshort_p,       [:M,                   ]),
    (:Cint,      :mpz_odd_p,               [:M,                   ]),
    (:Cint,      :mpz_even_p,              [:M,                   ]),
    (:Csize_t,   :mpz_sizeinbase,          [:M,   :I              ]),
                                                                   
    # 5.16 Special Functions                                       
    (:Csize_t,   :mpz_size,                [:M                    ]),
    (:mp_limb_t, :mpz_getlimbn,            [:M,   :MS             ]),
]



# --- Generation Loop ---
function gmp_ffi_defining_functions()
    for (ret_type, fname, arg_keys) in GMP_DEFINITIONS
        
        # シンボルが存在しない（純粋なマクロなど）場合はスキップ
        if !has_gmp_symbol(fname)
            @warn "GMP symbol __g$fname not found. Skipping auto-definition."
            continue
        end

        j_types = [TYPE_MAP[k][1] for k in arg_keys]
        c_types = [TYPE_MAP[k][2] for k in arg_keys]
        args    = [Symbol(:a, i)  for i in 1:lastindex(arg_keys)]
        
        # a1::T1, a2::T2...
        sig = [Expr(:(::), args[i], j_types[i]) for i in 1:lastindex(args)]
        
        # Cvoid の場合は Julia 側で Nothing として注釈を付ける
        j_ret = (ret_type === :Cvoid) ? :Nothing : ret_type
        
        gmp_target = (Symbol(:__g, fname), libgmp)

        @eval begin
            # 戻り値の型 ::$j_ret を追加
            export $fname
            function $fname($(sig...))::$j_ret
                return ccall($gmp_target, $ret_type, ($(c_types...),), $(args...))
            end
        end
    end
end

"""
    has_gmp_symbol(fname::Symbol)

libgmp 内に `__gmpz_` + fname というシンボルが存在するかを確認する。
"""
function has_gmp_symbol(fname::Symbol)
    # libgmp のハンドルを取得（すでにロードされているはずだが安全のために dlopen）
    lib_handle = dlopen(libgmp)
    
    # 実際のシンボル名（__gmpz_...）を生成
    gmp_sym_name = Symbol(:__g, fname)
    
    # シンボルを探す。見つかればそのアドレス（Ptr）、なければ C_NULL を返す
    ptr = dlsym_e(lib_handle, gmp_sym_name)
    
    return ptr != C_NULL
end

gmp_ffi_defining_functions()

export mpz_mod_ui, mpz_kronecker
mpz_mod_ui(a1::mpz_t, a2::mpz_t, a3::CulongMax)::Culong = mpz_fdiv_r_ui(a1, a2, a3)
mpz_kronecker(a1::mpz_t, a2::mpz_t)::Cint = mpz_jacobi(a1, a2)
