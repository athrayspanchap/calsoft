; IDL procedure by Athiray
; Copyright (c) 2017, FOXSI Mission University of Minnesota.  All rights reserved.
;       Unauthorized reproduction is allowed.


; Start		: 21 Jul 2017 15:56
; Last Mod 	: 31 Oct 2017 16:35

;-------------------  Details of the program --------------------------;
function gauss2,x,par
gauss_fit =(par[0]*exp(-0.5*((x-par[1])/par[2])^2))+(par[3]*exp(-0.5*((x-par[4])/par[5])^2))
return,gauss_fit
end



pro fwhm_gauss,el

if (el eq 'cr') then begin
    files='spec_data_dt6_cr.dat'
    restore,files
    peaks=[5.0]
    spec=spec_cr
endif

if (el eq 'fe') then begin
    files='spec_data_dt6_fe.dat'
    restore,files
    peaks=5.9
    spec=spec_fe
endif

if (el eq 'ni') then begin
    files='spec_data_dt6_ni.dat'
    restore,files
    peaks=7.4
    spec=spec_ni
endif
if (el eq 'cu') then begin
    files='spec_data_dt6_cu.dat'
    restore,files
    peaks=8.0
    spec=spec_cu
endif

if (el eq 'am') then begin
    files='spec_data_dt6_am.dat'
    restore,files
    ;peaks=[9.4,11.0,13.9,17.6,20.6,26.3,59.6]
    peaks=[13.9,17.7,26.3,59.6]
    ;peaks=[13.9,26.3,59.6]
    spec=spec_am
endif


;restore,'spec_data_dt6_ni.dat'
;spec=spec_ni

for asic=2,3 do  begin
    fname='efwhm_cal_asic'+strtrim(string(asic),2)+'.txt'
    openw,1,fname,/append

for k=0,n_elements(peaks)-1 do begin
    peak=peaks[k]
    ;peak=[7.46]
    bestsig=fltarr(100)
    no_iter=sqrt(n_elements(bestsig))
    sel_fwhm=fltarr(no_iter,no_iter)
    sel_peak=fltarr(no_iter,no_iter)
    step=0.04
    basic_range=[peak-0.1,peak+0.2]
    for i = 0,no_iter-1 do begin
	for j =0,no_iter-1 do begin
	    sel_range=[(basic_range[0]-(step*i)),(basic_range[1]+(step*j))]
	    print,sel_range
	    ;sel_fwhm[i,j]=calcfwhm(spec, asic=asic, ch=64, range=sel_range, rebin=1)
	    val=calcfwhm_psa(spec, asic=asic, ch=64, range=sel_range, rebin=1)
	    sel_fwhm[i,j]=val[1]
	    sel_peak[i,j]=val[0]
	    ;sel_fwhm[i,j]=fwhm
	    ;stop
	endfor
    endfor
    nan_in = where(finite(sel_fwhm) eq 0)
    sel_fwhm[nan_in] =1e5
    wr=where(sel_fwhm lt 0.1)
    sel_fwhm[wr]=1e5
    g_fwhm=min(sel_fwhm)
    in=where(sel_fwhm eq g_fwhm)
    g_peak=sel_peak[in]
    sig=g_fwhm/2.35
    peak=g_peak[0]
    en=spec[*,64,asic,0]
    sp=spec[*,64,asic,1]
;peaks_fe_ni_cu_am_cr_gauss,histfe=h_fe,histamni=h_ni,histamcu=h_cu,histame=h_am,histcr=h_cr
    sel_index=where((en ge (peak-(2*sig))) and (en le (peak+(2*sig))))
    ;sel_index=where((en ge (peak-(2.5*sig))) and (en le (peak+(2.5*sig))))
    selen=en[sel_index]
    selsp=sp[sel_index]
    pvalin=where(abs(selen-peak) eq min(abs(selen-peak)),numb)
    if numb eq 0 then stop
    pval=selsp[pvalin[0]]



    est=[pval,peak,sig]
    par=est
    p = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},(3))
    ;p(1).fixed=1
    ;p[*]=par
    if ((el eq 'fe') or (el eq 'cr')) then limit_fac=0.01 else limit_fac=0.08
    if ((el eq 'ni')) then limit_Fac=0.03
    if (peaks[k] eq 17.7) then limit_fac=0.03
    ;values for det. post 6
    if ((el eq 'fe') or (el eq 'cr')) then limit_fac=0.02 else limit_fac=0.08
    if ((el eq 'ni')) then limit_Fac=0.03
    if (peaks[k] eq 17.7) then limit_fac=0.03

    p(0).limited(0)=1
    p(0).limits(0)=par[0]-sqrt(par[0])
    p(0).limited(1)=1
    p(0).limits(1)=par[0]+sqrt(par[0])
    p(2).limited(0)=1
    p(2).limits(0)=par[2]-0.02
    p(2).limited(1)=1
    p(2).limits(1)=par[2]+limit_Fac



    ;fit=mpfitpeak(selen,selsp,coef,estimates=est,parinfo=p,nterms=3,sigma=si,measure_errors=sqrt(selsp))
    fit=mpfitpeak(selen,selsp,coef,estimates=est,parinfo=p,nterms=3,sigma=si,measure_errors=(0.1*selsp),chisq=chi)
    plot,selen,selsp,psym=10
    oplot,selen,fit
    print,coef,si
    another_gauss='n'
    ;READ, another_gauss, PROMPT=' Do you want to fit another gauss (y or n)? '

    if(another_gauss eq 'y') then begin
	peak2=''
	READ, peak2, PROMPT='enter peak energy in kev? '
	selen2_low=peak2-0.15
	selen2_high=peak2+0.3
	sel_index2=where((en ge selen2_low) and (en le selen2_high))
	selen2=en[sel_index2]
	selsp2=sp[sel_index2]
	est2=[1,peak2,sig]
	par2=est2
	p2 = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},(3))
	;p2(1).fixed=1
	;p2[1]=par2[1]
	p2(2).limited(0)=1
	p2(2).limits(0)=par2[2]-0.02
	p2(2).limited(1)=1
	p2(2).limits(1)=par2[2]+0.1




	fit2=mpfitpeak(selen2,selsp2,coef2,estimates=est2,nterms=3,sigma=si2,measure_errors=sqrt(selsp2),parinfo=p2)
	gpar=[coef,coef2]

        gparinf=replicate({value:0.D,fixed:0, limited:[0,0], limits:[0.D,0.D]},(3*2))
	gparinf(2).limited(0)=1
	gparinf(2).limits(0)=coef[2]-0.02
	gparinf(2).limited(1)=1
	gparinf(2).limits(1)=coef[2]+0.0
	gparinf(5).limited(0)=1
	gparinf(5).limits(0)=coef2[2]-0.02
	gparinf(5).limited(1)=1
	gparinf(5).limits(1)=coef2[2]+0.02

	gparinf(1).fixed= 1
    	gparinf[*].value=gpar


	;gparinf=[p,p2]
	en_low=min(selen)
	en_high=max(selen2)
	sel_index3=where((en ge en_low) and (en le en_high))
	en3=en[sel_index3]
	sp3=sp[sel_index3]
	er3=sqrt(sp3)
        gresult=mpfitfun('gauss2',en3,sp3,er3,gpar,perror=sig3,parinfo=gparinf);,sigma)
	ymodel = $
	(gresult[0]*exp(-0.5*((en3-gresult[1])/gresult[2])^2))+(gresult[3]*exp(-0.5*((en3-gresult[4])/gresult[5])^2))
   	 g1=(gresult[0]*exp(-0.5*((en3-gresult[1])/gresult[2])^2))
   	 g2=(gresult[3]*exp(-0.5*((en3-gresult[4])/gresult[5])^2))


    endif
    print,peak,g_fwhm,coef[1],coef[2]*2.35
    printf,1,coef[1],coef[2]*2.35,si[2]*2.35,chi,min(selen),max(selen)
   ;stop
endfor
close,1
endfor
;;print





;+
; NAME:
;      Function
;
; PURPOSE:
;       This function writes
;
; CATEGORY:
;       Input/Output
;
; CALLING SEQUENCE:
;       Function, input,output
END

