*----------------------------------------------------------------------
* sim10_3grp.do
*----------------------------------------------------------------------

clear all
version 14

capture mkdir "E:\ArielStuff\wtd_k_simulations"

local simreps   2000
local testreps  1000
local ess_ratio 0.9

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


capture program drop simreps_tail
program define simreps_tail, rclass
	args n0 n1 n2 ess_ratio testreps

	clear
	local ntot = `n0' + `n1' + `n2'
	quietly set obs `ntot'

	quietly gen byte group = 0            if _n <= `n0'
	quietly replace group  = 1            if _n >  `n0' & _n <= `n0'+`n1'
	quietly replace group  = 2            if _n >  `n0'+`n1'
	quietly gen double p   = runiform()

	mata: st_store(., st_addvar("double","gapval"), sim_gapfun(st_data(.,"p"), 50, 0.05, 0.95, 0.02))
	quietly gen double y = 50 + 10*invnormal(p) + (group==2)*gapval

	mata: st_store(., st_addvar("double","wgt"), sim_make_weights(`ntot', `ess_ratio'))

	quietly kstestk y [pweight=wgt], by(group) reps(`testreps')
	return scalar p_ks = r(p)

	quietly adtestk y [pweight=wgt], by(group) reps(`testreps')
	return scalar p_ad = r(p)

	quietly cvmtestk y [pweight=wgt], by(group) reps(`testreps')
	return scalar p_cvm = r(p)
end

capture program drop simreps_diffuse
program define simreps_diffuse, rclass
	args n0 n1 n2 ess_ratio testreps

	clear
	local ntot = `n0' + `n1' + `n2'
	quietly set obs `ntot'

	quietly gen byte group = 0            if _n <= `n0'
	quietly replace group  = 1            if _n >  `n0' & _n <= `n0'+`n1'
	quietly replace group  = 2            if _n >  `n0'+`n1'
	quietly gen double p   = runiform()

	quietly gen double y = 50 + 10*invnormal(p) + (group==2)*1

	mata: st_store(., st_addvar("double","wgt"), sim_make_weights(`ntot', `ess_ratio'))

	quietly kstestk y [pweight=wgt], by(group) reps(`testreps')
	return scalar p_ks = r(p)

	quietly adtestk y [pweight=wgt], by(group) reps(`testreps')
	return scalar p_ad = r(p)

	quietly cvmtestk y [pweight=wgt], by(group) reps(`testreps')
	return scalar p_cvm = r(p)
end


tempfile results_sens
local first = 1

local condlist tail diffuse
local ncond : word count `condlist'

forvalues c = 1/`ncond' {
	local cond : word `c' of `condlist'

	di as txt _n "=== re=0.9, `cond' condition: n0=600 n1=600 n2=800 (total=2000) ==="

	simulate p_ks=r(p_ks) p_ad=r(p_ad) p_cvm=r(p_cvm), ///
		reps(`simreps') seed(195`c'): ///
		simreps_`cond' 600 600 800 `ess_ratio' `testreps'

	quietly gen byte rej_ks  = (p_ks  < 0.05)
	quietly gen byte rej_ad  = (p_ad  < 0.05)
	quietly gen byte rej_cvm = (p_cvm < 0.05)
	quietly gen long  ntot   = 2000
	quietly gen double n0    = 600
	quietly gen double n1    = 600
	quietly gen double n2    = 800
	quietly gen double re    = 0.9
	quietly gen str12 condition = "`cond'"

	quietly summarize rej_ks, meanonly
	local rate_ks = r(mean)
	quietly summarize rej_ad, meanonly
	local rate_ad = r(mean)
	quietly summarize rej_cvm, meanonly
	local rate_cvm = r(mean)

	di as txt "  kstestk (KS)      rejection rate = " as res %5.3f `rate_ks'
	di as txt "  adtestk (AD)      rejection rate = " as res %5.3f `rate_ad'
	di as txt "  cvmtestk (CVM)    rejection rate = " as res %5.3f `rate_cvm'

	quietly keep condition re ntot n0 n1 n2 rej_ks rej_ad rej_cvm p_ks p_ad p_cvm
	if `first' {
		save `results_sens', replace
		local first = 0
	}
	else {
		append using `results_sens'
		save `results_sens', replace
	}
}

use `results_sens', clear
save "E:\ArielStuff\wtd_k_simulations\sim10_3grp.dta", replace
di as txt _n "Saved: E:\ArielStuff\wtd_k_simulations\sim10_3grp.dta"
