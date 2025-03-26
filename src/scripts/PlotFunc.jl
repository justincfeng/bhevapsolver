
#-----------------------------------------------------------------------
#   Scientific Notation Strings
#-----------------------------------------------------------------------
function sns( r , l=3 )   # General function
    n1=string(round(Float64(r),sigdigits=l))[1:l+1]
    n2=string(Int(floor(log10(r))))
    return "$n1 × 10^$n2"
end     #---------------------------------------------------------------

function snsl( r , l=3 , pref="" , suff="" )   # LaTeXStrings version
    n1=string(round(Float64(r),sigdigits=l))[1:l+1]
    n2=string(Int(floor(log10(r))))
    s1=pref*n1
    return L"%$s1 \times 10^{%$n2}%$suff"
end     #---------------------------------------------------------------

function snsM( r , l=3 )   # Mass label
    n1=string(round(Float64(r),sigdigits=l))[1:l+1]
    n2=string(Int(floor(log10(r))))
    return L"M/(%$n1 \times 10^{%$n2} M_\odot)"
end     #---------------------------------------------------------------

function snsMu( r , l=3 )   # Mass label
    n1=string(round(Float64(r),sigdigits=l))[1:l+1]
    n2=string(Int(floor(log10(r))))
    return L"%$n1 \times 10^{%$n2} M_\odot"
end     #---------------------------------------------------------------

function snsMh( r , l=3 )   # Mass label
    n1=string(round(Float64(r),sigdigits=l))[1:l+1]
    n2=string(Int(floor(log10(r))))
    return L"M_{\mathrm{h}}/(%$n1 \times 10^{%$n2} M_\odot)"
end     #---------------------------------------------------------------

function snst( r , l=3 )   # Time label
    n1=string(round(Float64(r),sigdigits=l))[1:l+1]
    n2=string(Int(floor(log10(r))))
    return L"t/(%$n1 \times 10^{%$n2} \mathrm{yr})"
end     #---------------------------------------------------------------

function snstv( r , l=3 )   # Time label
    n1=string(round(Float64(r),sigdigits=l))[1:l+1]
    n2=string(Int(floor(log10(r))))
    return L"t_{\mathrm{evap}}/(%$n1 \times 10^{%$n2} \mathrm{yr})"
end     #---------------------------------------------------------------
