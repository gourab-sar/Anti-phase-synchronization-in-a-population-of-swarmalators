!!!! determination of the states::

implicit none
integer,parameter::n=100,m=100,nt=50000,tr=45000
complex, parameter:: r= (0,1)   !sqrt(-1)  
real,parameter:: pi=acos(-1.0),pi2=2.0*pi
integer::i,j,i1,mm,i2,cpu_time
complex::c1,c2,dm1,dm2,dp1,dp2,an1,an2
real::theta0(1000),theta(1000),x0(1000),y0(1000),x(1000),y(1000),kt1(n),kt2(n),kt3(n),kt4(n),&
kx1(n),kx2(n),kx3(n),kx4(n),ky1(n),ky2(n),ky3(n),ky4(n),phi0(1000),phi(1000),zx0(1000),zy0(1000),&
zx(1000),zy(1000),kzx1(m),kzx2(m),kzx3(m),kzx4(m),kzy1(m),kzy2(m),kzy3(m),kzy4(m),kp1(m),kp2(m),&
kp3(m),kp4(m),sumsin(1000),sum1(1000),sum2(1000),sum3(1000),sum4(1000),sumx(1000),sumy(1000),&
sumzx(1000),sumzy(1000),sumt1(1000),sumt2(1000),sump1(1000),sump2(1000),a1,a2,b1,b2,a3,b3,&
distance,s1,s2,s3,zita(nt,1000),sumr1,sumr2,sumr3,sumr4,tem,gamm1,gamm2,gamm,tol,temp(nt),yo(nt),&
splus1,smin1,splus2,smin2,splus3,smin3,t,h,j1,j2,j3,k1,k2,k3


call srand(cpu_time)

h=0.01;

b1=1.0;b2=1.0;			!!!! repulsion co-eff
a1=1.0;a2=1.0;          !!!! attraction co-eff

!!!anti-phase sync parameters::
j1=0.1;k1=-0.1;   
j2=0.1;k2=-0.2; 

k3=-0.6;           !!!! phase coupling between two community
a3=1.0 ;            !!!! spatial attraction between two community
b3=1.0;			    !!!!repulsive intensity
j3=0.1;				!!!influence of phase on spatial coupling
tol=0.01;
temp=0.0;gamm1=0.0;gamm2=0.0;gamm=0.0;

open(10,file='x.dat')
open(20,file='y.dat')
open(30,file='theta.dat')
open(40,file='zx.dat')
open(50,file='zy.dat')
open(60,file='phi.dat')

!open(500,file='r1,-0.7.dat')
!open(600,file='r2,-0.7.dat')


!!!!___ initial condition____
do i=1,n
	x0(i)=2.0*rand()-1.0
	y0(i)=2.0*rand()-1.0
	theta0(i)=pi2*rand()
enddo
do i=1,m
	zx0(i)=2.0*rand()-1.0
	zy0(i)=2.0*rand()-1.0
	phi0(i)=pi2*rand() 
 enddo
 
 t=0.0;
 sumr1=0.0;sumr2=0.0;sumr3=0.0;sumr4=0.0;
 splus1=0.0; smin1=0.0
 splus2=0.0; smin2=0.0
 splus3=0.0; smin3=0.0
 
do mm=1,nt      !!!time loop
	print*,mm
	
!!!!1st
sumx=0.0;sumy=0.0;sumt1=0.0;sumt2=0.0;sum1=0.0;sum2=0.0;
sumzx=0.0;sumzy=0.0;sump1=0.0;sump2=0.0;sum3=0.0;sum4=0.0;
!!!___community n___
do i=1,n
	do j=1,n
		if(j.ne.i)then
			distance=sqrt((x0(j)-x0(i))**2.0+(y0(j)-y0(i))**2.0)
			sumt1(i)=sumt1(i)+sin(theta0(j)-theta0(i))/distance
			sumx(i)=sumx(i)+((x0(j)-x0(i))*(a1+j1*cos(theta0(j)-theta0(i))))/distance-b1*((x0(j)-x0(i))/distance**2.0)
			sumy(i)=sumy(i)+((y0(j)-y0(i))*(a1+j1*cos(theta0(j)-theta0(i))))/distance-b1*((y0(j)-y0(i))/distance**2.0)
		endif
	enddo
	
	do j=1,m
		distance=sqrt((zx0(j)-x0(i))**2.0+(zy0(j)-y0(i))**2.0)
		sumt2(i)=sumt2(i)+sin(phi0(j)-theta0(i))/distance
		sum1(i)=sum1(i)+((a3+j3*cos(phi0(j)-theta0(i)))*(zx0(j)-x0(i)))/distance-b3*((zx0(j)-x0(i))/distance**2.0)
		sum2(i)=sum2(i)+((a3+j3*cos(phi0(j)-theta0(i)))*(zy0(j)-y0(i)))/distance-b3*((zy0(j)-y0(i))/distance**2.0)
	enddo
	
		kx1(i)=sumx(i)/n+sum1(i)/m
		ky1(i)=sumy(i)/n+sum2(i)/m
		kt1(i)=(k1/n)*sumt1(i)+(k3/m)*sumt2(i)
enddo
   
!!!___community m___
do i=1,m
	do j=1,m
		if(j.ne.i)then
			distance=sqrt((zx0(j)-zx0(i))**2.0+(zy0(j)-zy0(i))**2.0)
			sump1(i)=sump1(i)+sin(phi0(j)-phi0(i))/distance
			sumzx(i)=sumzx(i)+((zx0(j)-zx0(i))*(a2+j2*cos(phi0(j)-phi0(i))))/distance-b2*((zx0(j)-zx0(i))/distance**2.0)
			sumzy(i)=sumzy(i)+((zy0(j)-zy0(i))*(a2+j2*cos(phi0(j)-phi0(i))))/distance-b2*((zy0(j)-zy0(i))/distance**2.0)
		endif
	enddo
	
	do j=1,n
		distance=sqrt((x0(j)-zx0(i))**2.0+(y0(j)-zy0(i))**2.0)
		sump2(i)=sump2(i)+sin(theta0(j)-phi0(i))/distance
		sum3(i)=sum3(i)+((a3+j3*cos(theta0(j)-phi0(i)))*(x0(j)-zx0(i)))/distance-b3*((x0(j)-zx0(i))/distance**2.0)
		sum4(i)=sum4(i)+((a3+j3*cos(theta0(j)-phi0(i)))*(y0(j)-zy0(i)))/distance-b3*((y0(j)-zy0(i))/distance**2.0)
	enddo
	
	kzx1(i)=sumzx(i)/m+sum3(i)/n
	kzy1(i)=sumzy(i)/m+sum4(i)/n
	kp1(i)=(k2/m)*sump1(i)+(k3/n)*sump2(i)
enddo

do i=1,n		
	x(i)=x0(i)+h*kx1(i)
	y(i)=y0(i)+h*ky1(i)
	theta(i)=theta0(i)+h*kt1(i)
enddo
do i=1,m
	zx(i)=zx0(i)+h*kzx1(i)
	zy(i)=zy0(i)+h*kzy1(i)
	phi(i)=phi0(i)+h*kp1(i)
enddo


!!!2nd
sumx=0.0;sumy=0.0;sumt1=0.0;sumt2=0.0;sum1=0.0;sum2=0.0;
sumzx=0.0;sumzy=0.0;sump1=0.0;sump2=0.0;sum3=0.0;sum4=0.0;
!!!!___community n___
do i=1,n
	do j=1,n
		if(j.ne.i)then
			distance=sqrt((x(j)-x(i))**2.0+(y(j)-y(i))**2.0)
			sumt1(i)=sumt1(i)+sin(theta(j)-theta(i))/distance
			sumx(i)=sumx(i)+((x(j)-x(i))*(a1+j1*cos(theta(j)-theta(i))))/distance-b1*((x(j)-x(i))/distance**2.0)
			sumy(i)=sumy(i)+((y(j)-y(i))*(a1+j1*cos(theta(j)-theta(i))))/distance-b1*((y(j)-y(i))/distance**2.0)
		endif
	enddo
	
	do j=1,m
		distance=sqrt((zx(j)-x(i))**2.0+(zy(j)-y(i))**2.0)
		sumt2(i)=sumt2(i)+sin(phi(j)-theta(i))/distance
		sum1(i)=sum1(i)+((a3+j3*cos(phi(j)-theta(i)))*(zx(j)-x(i)))/distance-b3*((zx(j)-x(i))/distance**2.0)
		sum2(i)=sum2(i)+((a3+j3*cos(phi(j)-theta(i)))*(zy(j)-y(i)))/distance-b3*((zy(j)-y(i))/distance**2.0)
	enddo
	
		kx2(i)=sumx(i)/n+sum1(i)/m
		ky2(i)=sumy(i)/n+sum2(i)/m
		kt2(i)=(k1/n)*sumt1(i)+(k3/m)*sumt2(i)
enddo
   
!!!___community m___
do i=1,m
	do j=1,m
		if(j.ne.i)then
			distance=sqrt((zx(j)-zx(i))**2.0+(zy(j)-zy(i))**2.0)
			sump1(i)=sump1(i)+sin(phi(j)-phi(i))/distance
			sumzx(i)=sumzx(i)+((zx(j)-zx(i))*(a2+j2*cos(phi(j)-phi(i))))/distance-b2*((zx(j)-zx(i))/distance**2.0)
			sumzy(i)=sumzy(i)+((zy(j)-zy(i))*(a2+j2*cos(phi(j)-phi(i))))/distance-b2*((zy(j)-zy(i))/distance**2.0)
		endif
	enddo
	
	do j=1,n
		distance=sqrt((x(j)-zx(i))**2.0+(y(j)-zy(i))**2.0)
		sump2(i)=sump2(i)+sin(theta(j)-phi(i))/distance
		sum3(i)=sum3(i)+((a3+j3*cos(theta(j)-phi(i)))*(x(j)-zx(i)))/distance-b3*((x(j)-zx(i))/distance**2.0)
		sum4(i)=sum4(i)+((a3+j3*cos(theta(j)-phi(i)))*(y(j)-zy(i)))/distance-b3*((y(j)-zy(i))/distance**2.0)
	enddo
	
	kzx2(i)=sumzx(i)/m+sum3(i)/n
	kzy2(i)=sumzy(i)/m+sum4(i)/n
	kp2(i)=(k2/m)*sump1(i)+(k3/n)*sump2(i)
enddo

!!!!total
do i=1,n		
	theta0(i)=theta0(i)+h*0.5*(kt1(i)+kt2(i))
	theta0(i)=modulo(theta0(i),pi2)
	x0(i)=x0(i)+h*0.5*(kx1(i)+kx2(i))
	y0(i)=y0(i)+h*0.5*(ky1(i)+ky2(i))
	zita(mm,i)=atan2(y0(i),x0(i))
enddo

do i=1,m
	phi0(i)=phi0(i)+h*0.5*(kp1(i)+kp2(i))
	phi0(i)=modulo(phi0(i),pi2)
	zx0(i)=zx0(i)+h*0.5*(kzx1(i)+kzx2(i))
	zy0(i)=zy0(i)+h*0.5*(kzy1(i)+kzy2(i))
	zita(mm,i)=atan2(zy0(i),zx0(i))
enddo

!!!_____order_parameter_____

 c1=0.0;c2=0.0;
 an1=0.0;an2=0.0;
 dm1=0.0;dm2=0.0;
 dp1=0.0;dp2=0.0;
 
if(mm.gt.tr)then 
  !!community_n
	do i=1,n
		 c1=c1+exp(r*theta0(i))
		 
		 dp1=dp1+exp(r*(atan2(y0(i),x0(i))+theta0(i)))
		 dm1=dm1+exp(r*(atan2(y0(i),x0(i))-theta0(i)))
		 an1=an1+exp(2*r*theta0(i))
		 !c3=c3+exp(r*theta0(i))
		 
	enddo
	
	!!community_m
	do i=1,m
		c2=c2+exp(r*phi0(i))
		
		dp2=dp2+exp(r*(atan2(zy0(i),zx0(i))+phi0(i)))
		dm2=dm2+exp(r*(atan2(zy0(i),zx0(i))-phi0(i)))
		an2=an2+exp(2*r*phi0(i))
		!c3=c3+exp(r*phi0(i))
	enddo
	
!!total order parameter
	sumr1=sumr1+cabs(c1)/n
	splus1=splus1+cabs(dp1)/n
	smin1=smin1+cabs(dm1)/n
	
	sumr2=sumr2+cabs(c2)/m
	splus2=splus2+cabs(dp2)/m
	smin2=smin2+cabs(dm2)/m
	
	sumr3=sumr3+cabs(c1+c2)/(n+m)
	sumr4=sumr4+cabs(an1+an2)/(n+m)
	smin3=smin3+cabs(dm1+dm2)/(n+m)
	splus3=splus3+cabs(dp1+dp2)/(n+m)
	
endif
	
	t=t+h
  write(10,*) (x0(i),i=1,n)
  write(20,*) (y0(i),i=1,n)
  write(30,*) (theta0(i),i=1,n)
  write(40,*) (zx0(i),i=1,m)
  write(50,*) (zy0(i),i=1,m)
  write(60,*) (phi0(i),i=1,m)


enddo   !!! time

!!!!____________calculation of gamma___________
		do i = 1,n
      do mm = 1,nt
      	yo(mm) = zita(mm,i)
      end do
      do mm = (nt/3)+1,nt
      	temp(mm-(nt/3)) = sin(yo(mm)) 
      end do
      tem = (maxval(temp) - minval(temp))/2.0
      if (tem .gt. (1.0-tol))then
      	gamm1 = gamm1 + 1.0
      end if
    end do
      
    do i = 1,m
      do mm = 1,nt
      	yo(mm) = zita(mm,i)
      end do
      do mm = (nt/3)+1,nt
      	temp(mm-(nt/3)) = sin(yo(mm)) 
      end do
      tem = (maxval(temp) - minval(temp))/2.0
      if (tem .gt. (1.0-tol))then
      	gamm2 = gamm2 + 1.0
      end if
    end do
      gamm=(gamm1+gamm2)/(n+m)
	  
write(*,*) k3, sumr3/(nt-tr), max(smin3,splus3)/(nt-tr), gamm, sumr4/(nt-tr)

end
