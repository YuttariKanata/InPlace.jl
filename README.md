# InPlace.jl

[![Build Status](https://github.com/YuttariKanata/InPlace.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/YuttariKanata/InPlace.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/YuttariKanata/InPlace.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/YuttariKanata/InPlace.jl)

InPlace.jl provides mutating arithmetic helpers for Julia's `BigInt` and
`BigFloat`.  It wraps GMP/MPFR operations behind Julia-style `!` functions so
that callers can reuse result objects in tight loops.

```julia
using InPlace

x = BigInt(100)
y = BigInt(50)
r = BigInt(0)

add!(r, x, y)       # r = 150
mul!(r, r, 10)      # r = 1500
div!(r, r, 7)       # r = 214, truncating division
rem!(r, 1500, 7)    # r = 2
```

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/YuttariKanata/InPlace.jl")
```

## Design

The first argument is usually the destination object and is mutated in-place:

```julia
result = BigInt(0)
add!(result, a, b)  # result := a + b
add!(result, result, c)
```

For integer arguments, InPlace.jl uses GMP's scalar entry points when the value
fits the platform C type (`Clong` or `Culong`), even if the Julia value is a
wider type such as `Int64` or `UInt64`.  Values that do not fit are converted to
`BigInt` and handled by the generic GMP path.

For `BigFloat`, operations use the precision of the destination and accept MPFR
rounding modes where applicable:

```julia
x = BigFloat(0; precision=256)
set!(x, "2")

r = BigFloat(0; precision=256)
sqrt!(r, x; rounding=InPlace.MPFRRoundNearest)
```

## Common BigInt Operations

```julia
set!(r, x)
swap!(x, y)

add!(r, x, y)
sub!(r, x, y)
mul!(r, x, y)
div!(r, x, y; rounding=:trunc)  # :trunc, :floor, or :ceil
rem!(r, x, y; rounding=:trunc)
mod!(r, x, y)

addmul!(r, x, y)     # r += x*y
submul!(r, x, y)     # r -= x*y
mul_2exp!(r, x, n)   # r = x * 2^n
div_2exp!(r, x, n)

pow!(r, x, n)
powm!(r, x, e, m)
powm_sec!(r, x, e, m)

sqrt!(r, x)
sqrtrem!(root, rem, x)
root!(r, x, n)
rootrem!(root, rem, x, n)

gcd!(r, x, y)
lcm!(r, x, y)
gcdext!(g, s, t, x, y)
invert!(r, x, m)
remove!(r, x, factor)

fac_ui!(r, n)
bin_ui!(r, n, k)
fib_ui!(r, n)
lucnum_ui!(r, n)

and!(r, x, y)
or!(r, x, y)
xor!(r, x, y)
com!(r, x)
setbit!(r, bit)
clrbit!(r, bit)
combit!(r, bit)
```

## Common BigFloat Operations

```julia
set!(r, x; rounding=InPlace.MPFRRoundNearest)
swap!(x, y)

add!(r, x, y; rounding=InPlace.MPFRRoundNearest)
sub!(r, x, y; rounding=InPlace.MPFRRoundNearest)
mul!(r, x, y; rounding=InPlace.MPFRRoundNearest)
div!(r, x, y; rounding=InPlace.MPFRRoundNearest)

sqrt!(r, x)
cbrt!(r, x)
rec_sqrt!(r, x)
square!(r, x)
pow!(r, x, y)
rootn!(r, x, n)

exp!(r, x)
exp2!(r, x)
exp10!(r, x)
log!(r, x)
log2!(r, x)
log10!(r, x)

sin!(r, x)
cos!(r, x)
tan!(r, x)
sin_cos!(s, c, x)

gamma!(r, x)
lngamma!(r, x)
zeta!(r, x)
erf!(r, x)
erfc!(r, x)

fma!(r, x, y, z)
fms!(r, x, y, z)
hypot!(r, x, y)
ceil!(r, x)
floor!(r, x)
round!(r, x)
trunc!(r, x)
frac!(r, x)
modf!(int_part, frac_part, x)

pi!(r)
euler!(r)
catalan!(r)
```

## Exported Functions

All public mutating helpers currently exported by InPlace.jl:

```julia
set!, swap!,
add!, sub!, mul!, div!, rem!, mod!,
neg!, abs!, min!, max!,

pow!, powm!, powm_sec!,
sqrt!, sqrtrem!, root!, rootrem!, rootn!,
square!, rec_sqrt!, cbrt!,

addmul!, submul!,
mul_2exp!, div_2exp!,
divexact!, divrem!, divrem_ui!,

gcd!, gcd_ui!, gcdext!, lcm!, lcm_ui!,
invert!, remove!, nextprime!,

fac_ui!, bin_ui!, bin_uiui!,
fib_ui!, fib2_ui!, lucnum_ui!, lucnum2_ui!,

and!, or!, xor!, com!,
setbit!, clrbit!, combit!,

sin!, cos!, tan!, sin_cos!,
asin!, acos!, atan!, atan2!,
sec!, csc!, cot!,

sinh!, cosh!, tanh!, sinh_cosh!,
asinh!, acosh!, atanh!,
sech!, csch!, coth!,

exp!, exp2!, exp10!, expm1!,
log!, log2!, log10!, log1p!,

gamma!, gamma_inc!, lngamma!, lgamma!, digamma!,
beta!, zeta!, erf!, erfc!,
j0!, j1!, jn!, y0!, y1!, yn!,
agm!, ai!, eint!, li2!,

rint!, ceil!, floor!, round!, trunc!,
frac!, modf!, fmodquo!, remquo!,
fma!, fms!, hypot!, reldiff!, dim!,

nexttoward!, copysign!, nextabove!, nextbelow!,
get_exp!, set_exp!, setsign!,

pi!, euler!, catalan!
```

## Notes

- Mutating functions return the destination object unless the underlying GMP or
  MPFR function naturally returns a status value.
- `invert!` throws an error when the modular inverse does not exist.
- `div!` and `rem!` on `BigInt` use `:trunc` by default; pass `:floor` or
  `:ceil` when that is the intended GMP division family.
- `BigFloat` conversions from decimal strings are preferable when exact decimal
  input matters. A literal such as `2.7` is already a `Float64`.

## License

InPlace.jl is licensed under the MIT License. See [LICENSE](LICENSE).

---

## Supported Low-level C APIs

This library provides Julia-side wrappers for the following functional groups defined in `gmp.h` and `mpfr.h`. Each Julia function `f!(rop, ...)` corresponds to its respective C counterpart (e.g., `mpz_f`, `mpfr_f`).

### GNU MP (GMP) Supported Functions

```julia
:mpz_realloc2        
:mpz_set             :mpz_set_ui          :mpz_set_si          :mpz_set_d           :mpz_set_str         
:mpz_swap            
:mpz_get_ui          :mpz_get_si          :mpz_get_d           :mpz_get_d_2exp      :mpz_get_str         
:mpz_add             :mpz_add_ui          :mpz_sub             :mpz_sub_ui          :mpz_ui_sub          
:mpz_mul             :mpz_mul_si          :mpz_mul_ui          :mpz_addmul          :mpz_addmul_ui       :mpz_submul          :mpz_submul_ui       :mpz_mul_2exp        
:mpz_neg             :mpz_abs             
:mpz_cdiv_q          :mpz_cdiv_r          :mpz_cdiv_qr         :mpz_cdiv_q_ui       :mpz_cdiv_r_ui       :mpz_cdiv_qr_ui      :mpz_cdiv_ui         :mpz_cdiv_q_2exp     :mpz_cdiv_r_2exp     
:mpz_fdiv_q          :mpz_fdiv_r          :mpz_fdiv_qr         :mpz_fdiv_q_ui       :mpz_fdiv_r_ui       :mpz_fdiv_qr_ui      :mpz_fdiv_ui         :mpz_fdiv_q_2exp     :mpz_fdiv_r_2exp     
:mpz_tdiv_q          :mpz_tdiv_r          :mpz_tdiv_qr         :mpz_tdiv_q_ui       :mpz_tdiv_r_ui       :mpz_tdiv_qr_ui      :mpz_tdiv_ui         :mpz_tdiv_q_2exp     :mpz_tdiv_r_2exp     
:mpz_mod             
:mpz_divexact        :mpz_divexact_ui     :mpz_divisible_p     :mpz_divisible_ui_p  :mpz_divisible_2exp_p
:mpz_congruent_p     :mpz_congruent_ui_p  :mpz_congruent_2exp_p
:mpz_powm            :mpz_powm_ui         :mpz_powm_sec        
:mpz_pow_ui          :mpz_ui_pow_ui       
:mpz_root            :mpz_rootrem         :mpz_sqrt            :mpz_sqrtrem         :mpz_perfect_power_p  :mpz_perfect_square_p 
:mpz_probab_prime_p  :mpz_nextprime       :mpz_prevprime       
:mpz_gcd             :mpz_gcd_ui          :mpz_gcdext          :mpz_lcm             :mpz_lcm_ui          
:mpz_invert          
:mpz_jacobi          :mpz_legendre        :mpz_kronecker_si    :mpz_kronecker_ui    :mpz_si_kronecker    :mpz_ui_kronecker    
:mpz_remove          
:mpz_fac_ui          :mpz_2fac_ui         :mpz_mfac_uiui       :mpz_primorial_ui    
:mpz_bin_ui          :mpz_bin_uiui        
:mpz_fib_ui          :mpz_fib2_ui         :mpz_lucnum_ui       :mpz_lucnum2_ui      
:mpz_cmp             :mpz_cmp_d           :mpz_cmp_si          :mpz_cmp_ui          :mpz_cmpabs          :mpz_cmpabs_d        :mpz_cmpabs_ui       
:mpz_and             :mpz_ior             :mpz_xor             :mpz_com             
:mpz_popcount        :mpz_hamdist         
:mpz_scan0           :mpz_scan1           
:mpz_setbit          :mpz_clrbit          :mpz_combit          :mpz_tstbit          
:mpz_out_str         :mpz_inp_str         :mpz_out_raw         :mpz_inp_raw         
:mpz_urandomb        :mpz_urandomm        :mpz_rrandomb        :mpz_random          :mpz_random2         
:mpz_import          :mpz_export          
:mpz_fits_ulong_p    :mpz_fits_slong_p    :mpz_fits_uint_p     :mpz_fits_sint_p     :mpz_fits_ushort_p   :mpz_fits_sshort_p   
:mpz_sizeinbase      :mpz_size            
:mpz_getlimbn        :mpz_limbs_read      :mpz_limbs_write     :mpz_limbs_modify    :mpz_limbs_finish    
```

#### GNU MP (GMP) Unsupported Functions

```julia
:mpz_init,           :mpz_init2,          :mpz_init_set,       :mpz_init_set_d,     :mpz_init_set_si,    :mpz_init_set_str,   :mpz_init_set_ui,    :mpz_inits,          :mpz_clear,          :mpz_clears,         :mpz_array_init,     :mpz_roinit_n,       :mpz_set_q,          :mpz_set_f
```

### GNU MPFR Supported Functions

```julia
:mpfr_set_default_prec             :mpfr_get_default_prec             :mpfr_set_default_rounding_mode    :mpfr_get_default_rounding_mode    
:mpfr_set_prec                     :mpfr_get_prec                     
:mpfr_prec_round                   :mpfr_can_round                    :mpfr_min_prec                     
:mpfr_set                          :mpfr_set_ui                       :mpfr_set_si                       :mpfr_set_uj                       :mpfr_set_sj                       :mpfr_set_flt                      :mpfr_set_d                        :mpfr_set_z                        :mpfr_set_ui_2exp                  :mpfr_set_si_2exp                  :mpfr_set_uj_2exp                  :mpfr_set_sj_2exp                  :mpfr_set_z_2exp                   :mpfr_set_str                      
:mpfr_strtofr                      
:mpfr_set_nan                      :mpfr_set_inf                      :mpfr_set_zero                     
:mpfr_swap                         
:mpfr_get_flt                      :mpfr_get_d                        :mpfr_get_si                       :mpfr_get_ui                       :mpfr_get_sj                       :mpfr_get_uj                       :mpfr_get_d_2exp                   
:mpfr_frexp                        
:mpfr_get_z_2exp                   :mpfr_get_z                        :mpfr_get_str_ndigits              :mpfr_get_str                      
:mpfr_free_str                     
:mpfr_fits_ulong_p                 :mpfr_fits_slong_p                 :mpfr_fits_uint_p                  :mpfr_fits_sint_p                  :mpfr_fits_ushort_p                :mpfr_fits_sshort_p                :mpfr_fits_uintmax_p               :mpfr_fits_intmax_p                
:mpfr_add                          :mpfr_add_ui                       :mpfr_add_si                       :mpfr_add_d                        :mpfr_add_z                        
:mpfr_sub                          :mpfr_ui_sub                       :mpfr_sub_ui                       :mpfr_si_sub                       :mpfr_sub_si                       :mpfr_d_sub                        :mpfr_sub_d                        :mpfr_z_sub                        :mpfr_sub_z                        
:mpfr_mul                          :mpfr_mul_ui                       :mpfr_mul_si                       :mpfr_mul_d                        :mpfr_mul_z                        
:mpfr_sqr                          
:mpfr_div                          :mpfr_ui_div                       :mpfr_div_ui                       :mpfr_si_div                       :mpfr_div_si                       :mpfr_d_div                        :mpfr_div_d                        :mpfr_div_z                        
:mpfr_sqrt                         :mpfr_sqrt_ui                      :mpfr_rec_sqrt                     :mpfr_cbrt                         :mpfr_rootn_ui                     :mpfr_rootn_si                     :mpfr_root                         
:mpfr_neg                          :mpfr_abs                          :mpfr_dim                          
:mpfr_mul_2ui                      :mpfr_mul_2si                      :mpfr_div_2ui                      :mpfr_div_2si                      :mpfr_mul_2exp                     :mpfr_div_2exp                     
:mpfr_fac_ui                       
:mpfr_fma                          :mpfr_fms                          :mpfr_fmma                         :mpfr_fmms                         
:mpfr_hypot                        
:mpfr_reldiff                      
:mpfr_min                          :mpfr_max                          
:mpfr_sum                          
:mpfr_dot                          
:mpfr_cmp                          :mpfr_cmp_ui                       :mpfr_cmp_si                       :mpfr_cmp_d                        :mpfr_cmp_z                        :mpfr_cmp_ui_2exp                  :mpfr_cmp_si_2exp                  :mpfr_cmpabs                       :mpfr_cmpabs_ui                    
:mpfr_nan_p                        :mpfr_inf_p                        :mpfr_number_p                     :mpfr_zero_p                       :mpfr_regular_p                    
:mpfr_sgn                          
:mpfr_greater_p                    :mpfr_greaterequal_p               :mpfr_less_p                       :mpfr_lessequal_p                  :mpfr_equal_p                      :mpfr_lessgreater_p                :mpfr_unordered_p                  :mpfr_total_order_p                
:mpfr_log                          :mpfr_log_ui                       :mpfr_log2                         :mpfr_log10                        :mpfr_log1p                        :mpfr_log2p1                       :mpfr_log10p1                      
:mpfr_exp                          :mpfr_exp2                         :mpfr_exp10                        :mpfr_expm1                        :mpfr_exp2m1                       :mpfr_exp10m1                      
:mpfr_pow                          :mpfr_powr                         :mpfr_pow_ui                       :mpfr_pow_si                       :mpfr_pow_uj                       :mpfr_pow_sj                       :mpfr_pown                         :mpfr_pow_z                        :mpfr_ui_pow_ui                    :mpfr_ui_pow                       :mpfr_compound_si                  
:mpfr_cos                          :mpfr_sin                          :mpfr_tan                          :mpfr_cosu                         :mpfr_sinu                         :mpfr_tanu                         :mpfr_cospi                        :mpfr_sinpi                        :mpfr_tanpi                        :mpfr_sin_cos                      :mpfr_sec                          :mpfr_csc                          :mpfr_cot                          :mpfr_acos                         :mpfr_asin                         :mpfr_atan                         :mpfr_acosu                        :mpfr_asinu                        :mpfr_atanu                        :mpfr_acospi                       :mpfr_asinpi                       :mpfr_atanpi                       :mpfr_atan2                        :mpfr_atan2u                       :mpfr_atan2pi                      :mpfr_cosh                         :mpfr_sinh                         :mpfr_tanh                         :mpfr_sinh_cosh                    :mpfr_sech                         :mpfr_csch                         :mpfr_coth                         :mpfr_acosh                        :mpfr_asinh                        :mpfr_atanh                        
:mpfr_eint                         
:mpfr_li2                          
:mpfr_gamma                        :mpfr_gamma_inc                    :mpfr_lngamma                      :mpfr_lgamma                       :mpfr_digamma                      
:mpfr_beta                         
:mpfr_zeta                         :mpfr_zeta_ui                      
:mpfr_erf                          :mpfr_erfc                         
:mpfr_j0                           :mpfr_j1                           :mpfr_jn                           :mpfr_y0                           :mpfr_y1                           :mpfr_yn                           
:mpfr_agm                          
:mpfr_ai                           
:mpfr_const_log2                   :mpfr_const_pi                     :mpfr_const_euler                  :mpfr_const_catalan                
:mpfr_rint                         :mpfr_ceil                         :mpfr_floor                        :mpfr_round                        :mpfr_roundeven                    :mpfr_trunc                        :mpfr_rint_ceil                    :mpfr_rint_floor                   :mpfr_rint_round                   :mpfr_rint_roundeven               :mpfr_rint_trunc                   
:mpfr_frac                         
:mpfr_modf                         :mpfr_fmod                         :mpfr_fmodquo                      :mpfr_remainder                    :mpfr_remquo                       
:mpfr_integer_p                    
:mpfr_nexttoward                   :mpfr_nextabove                    :mpfr_nextbelow                    
:mpfr_get_exp                      :mpfr_set_exp                      
:mpfr_signbit                      :mpfr_setsign                      :mpfr_copysign                     
:mpfr_out_str                      :mpfr_inp_str                      
:mpfr_fpif_export                  :mpfr_fpif_import                  
:mpfr_dump                         :mpfr_urandomb                     :mpfr_urandom                      :mpfr_nrandom                      :mpfr_erandom                      :mpfr_grandom                      
:mpfr_clear_flags                  :mpfr_clear_underflow              :mpfr_clear_overflow               :mpfr_clear_divby0                 :mpfr_clear_nanflag                :mpfr_clear_inexflag               :mpfr_clear_erangeflag             
:mpfr_underflow_p                  :mpfr_overflow_p                   :mpfr_divby0_p                     :mpfr_nanflag_p                    :mpfr_inexflag_p                   :mpfr_erangeflag_p                 
:mpfr_flags_save                   :mpfr_flags_restore                :mpfr_flags_clear                  :mpfr_flags_set                    :mpfr_flags_test                   
:mpfr_set_underflow                :mpfr_set_overflow                 :mpfr_set_divby0                   :mpfr_set_nanflag                  :mpfr_set_inexflag                 :mpfr_set_erangeflag               :mpfr_get_emin                     :mpfr_get_emax                     :mpfr_set_emin                     :mpfr_set_emax                     :mpfr_get_emin_min                 :mpfr_get_emin_max                 :mpfr_get_emax_min                 :mpfr_get_emax_max                 
:mpfr_check_range                  
:mpfr_subnormalize                 
:mpfr_get_version                  :mpfr_get_patches                  
:mpfr_buildopt_tls_p               :mpfr_buildopt_float128_p          :mpfr_buildopt_decimal_p           :mpfr_buildopt_gmpinternals_p      :mpfr_buildopt_sharedcache_p       :mpfr_buildopt_tune_case           
:mpfr_free_cache                   :mpfr_free_cache2                  :mpfr_free_pool                    
:mpfr_mp_memory_cleanup            
:mpfr_custom_get_size              :mpfr_custom_init                  :mpfr_custom_init_set              :mpfr_custom_get_significand       :mpfr_custom_get_exp               :mpfr_custom_get_kind              :mpfr_custom_move                  
```

#### GNU MPFR Unsupported Functions

```julia
:mpfr_get_ld                   :mpfr_set_f                    :mpfr_div_q                    :mpfr_init2                    :mpfr_set_decimal128           :mpfr_init_set_ld              :mpfr_init_set_ui              :mpfr_get_decimal64            :mpfr_init_set_d               :mpfr_cmp_q                    :mpfr_get_q                    :mpfr_get_decimal128           :mpfr_init_set                 :mpfr_cmp_f                    :mpfr_clear                    :mpfr_init                     :mpfr_set_float128             :mpfr_vsprintf                 :mpfr_inits                    :mpfr_add_q                    :mpfr_set_ld                   :mpfr_vasprintf                :mpfr_asprintf                 :mpfr_clears                   :mpfr_get_f                    :mpfr_fprintf                  :mpfr_init_set_q               :mpfr_init_set_z               :mpfr_mul_q                    :mpfr_vsnprintf                :mpfr_set_q                    :mpfr_set_decimal64            :mpfr_sub_q                    :mpfr_printf                   :mpfr_set_prec_raw             :mpfr_vfprintf                 :mpfr_inits2                   :mpfr_snprintf                 :mpfr_vprintf                  :mpfr_init_set_si              :mpfr_init_set_str             :mpfr_init_set_f               :mpfr_get_float128             :mpfr_get_ld_2exp              :mpfr_sprintf                  
```
