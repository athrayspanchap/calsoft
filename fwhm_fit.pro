; IDL procedure by Athiray
; Copyright (c) 2017, FOXSI Mission University of Minnesota.  All rights reserved.
;       Unauthorized reproduction is allowed.


; Start		: 23 Jul 2017 03:06
; Last Mod 	: 03 Aug 2017 11:34

;-------------------  Details of the program --------------------------;
function efwhmfit,x,p
w=0.0036
Fano=0.115
Ymod = sqrt(p[0]^2+(2.35^2*Fano*w*x)+(p[1]*x)^2)
return, ymod
end


PRO fwhm_fit

for asic=2,3 do begin
    fname='efwhm_cal_asic'+strtrim(string(asic),2)+'.txt'
    readcol,fname,energy,fwhm,fwhmerr
    fwhmerr = fwhmerr ;0.01*fwhm
    o=where(fwhmerr eq 0.)
    fwhmerr[o]=0.02*fwhm[o]
    para=[1.D,0.005D]
    initguess = para
    result=mpfitfun('efwhmfit',energy,fwhm,fwhmerr,para,perror=sig);,sigma)
    Fano=0.115
    w=0.0036
    ymodel = sqrt(result[0]^2+(2.35^2*Fano*w*energy)+(result[1]*energy)^2)
    ymin=    sqrt((result[0]-(sig[0]))^2+(2.35^2*Fano*w*energy)+((result[1]-(sig[1]))*energy)^2)
    ymax=    sqrt((result[0]+(sig[0]))^2+(2.35^2*Fano*w*energy)+((result[1]+(sig[1]))*energy)^2)
    pfname='efwhm_cal_asic'+strtrim(string(asic),2)+'.eps'
    set_plot,'PS'
    loadct,11
    fname='efwhm_cal_asic'+strtrim(string(asic),2)+'.eps'
    device,filename=fname,/color
    ploterr,energy,fwhm,fwhmerr,xtitle='Energy (keV)',ytitle='FWHM (keV)',$
	psym=1,thick=5,xthick=4,ythick=4,charsize=1.75,charthick=4,font=1
    oplot,energy(sort(energy)),ymodel(sort(energy)),color=250,thick=3
    oplot,energy(sort(energy)),ymin(sort(energy)),color=150,thick=3,linestyle=2
    oplot,energy(sort(energy)),ymax(sort(energy)),color=150,thick=3,linestyle=2
    asicnoise=string(result[0],format='(F6.4)')
    erasicnoise=string(sig[0],format='(F6.4)')
    edepnoise=string(result[1],format='(F6.4)')
    eredepnoise=string(sig[1],format='(F6.4)')
    cgtext,0.3,0.8,'ASIC '+strtrim(string(asic),2),/normal,color=120,$
		 font=1,orientation=0,charsize=2.52,charthick=3
    ;cgtext,0.2,0.75,edepnoise,/normal,$

	;font=1,orientation=0,charsize=1.2,charthick=3
    device,/close
    set_plot,'x'

    stop
endfor





END

