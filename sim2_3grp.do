*----------------------------------------------------------------------
* sim2_3grp.do
*----------------------------------------------------------------------

clear all
version 14

capture mkdir "E:\ArielStuff\wtd_k_simulations"

 local simreps   2000
 local testreps  1000

local ess_ratio 0.6

mata:

real vector sim_make_weights(real scalar n, real scalar ess_ratio)
{
	real scalar cv2, sigma2, mu
	real vector w

	if (ess_ratio >= 0.999) {
		return(J(n,1,1))
	}
	cv2    = 1/ess_ratio - 1
	sigma2 = ln(1+cv2)
	mu     = -sigma2/2
	w      = exp(rnormal(n,1,mu,sqrt(sigma2)))
	return(w)
}

real vector sim_gapfun(real vector p, real scalar gsize, real scalar p1, real scalar p2, real scalar width)
{
	real vector pa, ga, out
	real scalar i, j, n, m

	pa = (0.0005, p1-width, p1, p1+width, 0.50, p2-width, p2, p2+width, 0.9995)'
	ga = (0, 0, -gsize, 0, 0, 0, gsize, 0, 0)'

	n = rows(p)
	m = rows(pa)
	out = J(n, 1, .)
	for (i=1; i<=n; i++) {
		if (p[i] <= pa[1]) {
			out[i] = ga[1]
		}
		else if (p[i] >= pa[m]) {
			out[i] = ga[m]
		}
		else {
			for (j=1; j<m; j++) {
				if (p[i] >= pa[j] & p[i] <= pa[j+1]) {
					out[i] = ga[j] + (ga[j+1]-ga[j]) * (p[i]-pa[j])/(pa[j+1]-pa[j])
					j = m
				}
			}
		}
	}
	return(out)
}
end


capture program drop simreps_null
program define simreps_null, rclass
	args n0 n1 n2 ess_ratio testreps

	clear
	local ntot = `n0' + `n1' + `n2'
	quietly set obs `ntot'

	quietly gen byte group = 0            if _n <= `n0'
	quietly replace group  = 1            if _n >  `n0' & _n <= `n0'+`n1'
	quietly replace group  = 2            if _n >  `n0'+`n1'
	quietly gen double p   = runiform()

	quietly gen double y = 50 + 10*invnormal(p)

	mata: st_store(., st_addvar("double","wgt"), sim_make_weights(`ntot', `ess_ratio'))

	quietly kstestk y [pweight=wgt], by(group) reps(`testreps')
	return scalar p_ks = r(p)

	quietly adtestk y [pweight=wgt], by(group) reps(`testreps')
	return scalar p_ad = r(p)

	quietly cvmtestk y [pweight=wgt], by(group) reps(`testreps')
	return scalar p_cvm = r(p)
end

tempfile results8
local first = 1

local sizelist_n0 900 1200
local sizelist_n1 900 1200
local sizelist_n2 1200 1600
local nsizes : word count `sizelist_n0'

forvalues s = 1/`nsizes' {
	local n0  : word `s' of `sizelist_n0'
	local n1  : word `s' of `sizelist_n1'
	local n2  : word `s' of `sizelist_n2'
	local ntot = `n0' + `n1' + `n2'

	di as txt _n "=== NULL (Type I error / size check): n0=`n0' n1=`n1' n2=`n2' (total=`ntot') ==="

	simulate p_ks=r(p_ks) p_ad=r(p_ad) p_cvm=r(p_cvm), ///
		reps(`simreps') seed(215`s'): ///
		simreps_null `n0' `n1' `n2' `ess_ratio' `testreps'

	quietly gen byte rej_ks  = (p_ks  < 0.05)
	quietly gen byte rej_ad  = (p_ad  < 0.05)
	quietly gen byte rej_cvm = (p_cvm < 0.05)
	quietly gen long  ntot   = `ntot'
	quietly gen double n0    = `n0'
	quietly gen double n1    = `n1'
	quietly gen double n2    = `n2'
	quietly gen str12 condition = "null"

	quietly summarize rej_ks, meanonly
	local rate_ks = r(mean)
	quietly summarize rej_ad, meanonly
	local rate_ad = r(mean)
	quietly summarize rej_cvm, meanonly
	local rate_cvm = r(mean)

	di as txt "  kstestk (KS)      rejection rate = " as res %5.3f `rate_ks'
	di as txt "  adtestk (AD)      rejection rate = " as res %5.3f `rate_ad'
	di as txt "  cvmtestk (CVM)    rejection rate = " as res %5.3f `rate_cvm'

	quietly keep condition ntot n0 n1 n2 rej_ks rej_ad rej_cvm p_ks p_ad p_cvm
	if `first' {
		save `results8', replace
		local first = 0
	}
	else {
		append using `results8'
		save `results8', replace
	}
}

use `results8', clear
save "E:\ArielStuff\wtd_k_simulations\sim2_3grp.dta", replace
di as txt _n "Saved: E:\ArielStuff\wtd_k_simulations\sim2_3grp.dta"
