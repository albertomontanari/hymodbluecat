subroutine hymodfortran(param,area,tdelta,e,p,ntstep,qinitial,evapt,qtslow,qtquick,qtot)
!
!  Declarations for subroutine arguments.
!
implicit none
integer tdelta,ntstep,deltat,i,j
double precision er1,er2,e(ntstep),p(ntstep),qinitial,qtquick(ntstep),qtslow(ntstep), &
                 qtot(ntstep),evapt(ntstep),param(5),area,cmax,betap,alfa,kslow,kquick,w1,w2,fatconv, &
                 c1,c2,wquick(3),qquick,wslow,uslow,qslow,uquick,dummy

!  Parameters
	cmax=param(1)
	betap=param(2)
	alfa=param(3)
	kslow=param(4)
	kquick=param(5)

!  Wathershed area: conversion from kmq to mq
	area=area*1000000.

!  Set the conversion factor
	fatconv=(1./1000.0/tdelta)*area 
  
!  Store initialization
	w1=0.0 
	w2=0
	c1=0.0 
	c2=0.0 
	uquick=0.0 
	uslow=0.0 
 	wslow=qinitial/(1./kslow*fatconv)
	 do i=1,3
    	   wquick(i)=0.0
         end do

!---------------------------------------------------
!  Model loop  
!
	 do i=1,ntstep
!  Compute excess precipitation and evaporation
    	   w1=w2
           ! remove negative number
           	dummy=(1. - ((betap + 1.) * w1/cmax))
        	dummy=max(dummy,0.0)
           c1=cmax * (1. - (dummy**(1./(betap + 1.))))
	   c2=min(c1+p(i),cmax) 
	   c2=max(c2,0.0) 
	   er1=max((p(i) - cmax + c1),0.0)
	   	dummy=1. - c2/cmax
	   	dummy=max(dummy,0.0)
	   w2=(cmax/(betap + 1.))*(1. - (dummy**(betap + 1.)))
	   er2=max((c2 - c1) - (w2 - w1),0.0)
	   evapt(i)=(1.-(((cmax-c2)/(betap+ 1.))/(cmax/(betap+ 1.))))*e(i)
	   w2=max(w2 - evapt(i),0.0)  
    
!  Partition uquick and uslow into quick and slow flow component 
	   uquick=alfa*(er2+er1)
	   uslow=(1.- alfa)*(er2+er1)

!  Route slow flow component with single linear reservoir (kslow)
! The linear reservoir is resolved with the following relationship:
! w(t)=w(t-1)+p(t) - (w(t-1)+p(t))/k*Dt = w(t-1)+p(t) - (w(t-1)+p(t))/(k/Dt)
! If k is measured in units of Dt then Dt=1 and one obtains: 
! w(t) = (w(t-1)+p(t)) - (w(t-1)+p(t))/k = (w(t-i)+p(t))(1-1/k), then:
! w(t) = (1-1/k) w(t-1) + (1-1/k)p(t), which is the equation below.

	   wslow=(1. - 1./kslow)*wslow + (1.- 1./kslow)*uslow

! Then, from the above proof it follows that q(t) = (w(t-1)+p(t))/k
! Namely: 1/k/(1-1/k) w(t), which is the equation below.

	   qslow=(1./kslow/(1. - 1./kslow))*wslow
	   qslow=max(qslow,0.0) 
	        
! Quick flow # We repeat the above computation three times.

 	   qquick=0.0
	     do j=1,3
		wquick(j)=(1.- 1./kquick)*wquick(j) + (1. - 1./kquick)*uquick
		qquick=(1./kquick/(1.- 1./kquick))*wquick(j)
		uquick=qquick
 	     end do
 	   qquick=max(qquick,0.0)
	    
!  Compute quick, slow and total flow for timestep
 	   qtslow(i)=qslow*fatconv
	   qtquick(i)=qquick*fatconv
  	   qtot(i)=qtquick(i) + qtslow(i)
  
	 end do
!---------------------------------------------------
      end subroutine hymodfortran

