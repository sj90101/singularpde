  function inireg3j
% *****************************************************************
% * Produces top row in Figure 8 in                               *
% * "On integral equation methods for the first Dirichlet problem *
% *  of the biharmonic and modified biharmonic equations"         *
% * SIAM J. Sci. Comput., vol. 40(4), pp. A2609-A2630, 2018.      *
% * https://doi.org/10.1137/17M1162238                            *
% *****************************************************************  
  close all
  format longE
  format compact
% *** Comments on the parameters ************************************
%
% asymp: set asymp=1 if asymptotics is of prime concern
%        set asymp=0 if field evaluation is of prime concern
%
% iref:  Coarse mesh refinement towards corner vertex. Recommended
%        value is iref=0. A larger iref, in combination with dlim1=1,
%        could lead to speedups in the post-processor.
%
% *** User specified quantities *************************************
  theta=1.0*pi/2;  % corner opening angle
  t=theta*2/pi;    % tested case: t \in [0.4,1.7]
  if t<0.5
    npan0=42       % number of un-refined coarse panels
  elseif t<0.7
    npan0=38
  elseif t<1
    npan0=34
  else
    npan0=28
  end  
  nsub=40         % number of subdivisions
  asymp=0;        % see explanation above
  iref=0;         % Coarse mesh refinement towards corner vertex
  dlim1=0.8       % closeness parameter for 16-point quadrature
  ngr=300;        % field evaluation at ngr^2 points
  xylim=[-0.1 1.1 -0.5353 0.5353];  % computational domain
% ****************************************************************
  Ng=16;
  T=Tinit16;
  W=Winit16;
  [IP,IPW]=IPinit(T,W);
  Pbc=Pbcinit(IP);
  PWbc=Pbcinit(IPW);
  Pbc2=blkdiag(Pbc,Pbc);
  PWbc2=blkdiag(PWbc,PWbc);
%
% *** Panels, points, and weights ***
  [sinter,sinterdiff,npan]=panelinit(iref,npan0);
  [z,~,nz,wzp,kappa,awzp,zinter,np,~]=zinit(theta,sinter,sinterdiff,T,W,npan);
  ntot=2*np;
  [xg,yg,zg,ngrtot]=targinit(xylim,ngr);
%
% *** The M_coa^{\circ} and K_coa^{\circ} matrices are set up ***
  Mcirc=MKmatinit(z,nz,wzp,kappa,awzp,np);
  starind=[np-2*Ng+1:np 1:2*Ng 2*np-2*Ng+1:2*np np+1:np+2*Ng];
  myind1=[np-2*Ng+1:np 2*np-2*Ng+1:2*np ];
  myind2=[1:2*Ng np+1:np+2*Ng];
  Mcirc(myind1,myind2)=zeros(4*Ng);
  Mcirc(myind2,myind1)=zeros(4*Ng);  
%
% *** Construct the right hand side ***
  [g,gz]=testfun(z);
%  figure(1)
%  plot(real(z),imag(z))
%  hold on
%  plot(real(zs),imag(zs),'*')  
%  axis equal
%  drawnow
  gn=real(conj(nz).*gz);
  rhs=2*[g;gn];  
%  rhs=ones(2*np,1);  
%  
% *** Initializer ***
  starL=[Ng+1:5*Ng 7*Ng+1:11*Ng];
  circL=[1:Ng 5*Ng+1:7*Ng 11*Ng+1:12*Ng];
  t1=cputime;
  R0=initializeR0(T,W,theta,Pbc2,PWbc2,starL,circL);
  t2=cputime-t1;
  disp(['The CPU time on the initializer is ',num2str(t2), ' seconds'])
  
  if asymp==1
    prema=2
  else
    prema=min(nsub-ceil(log2(sum(awzp(1:Ng))/min(abs(zg)))),nsub)
  end
  nse=nsub-prema+1;
  disp(['nse = ',num2str(nse)])
  t1=cputime;
  [R,Rstor,Kstor]=Rcomp(R0,theta,T,W,Pbc2,PWbc2,nsub,prema,nse,npan0, ...
			iref,starL,circL);
  t2=cputime-t1;
  disp(['The CPU time on the recursion is ',num2str(t2), ' seconds'])
 
  RG=speye(ntot);
  RG(starind,starind)=R;  
%
% *** Solving main linear system ***
  t1=cputime;
  [rhotilde,iter]=myGMRESR(Mcirc,RG,rhs,ntot,100,2*eps);
  t2=cputime-t1;
  disp(['The CPU time on GMRES is ',num2str(t2), ' seconds'])

  disp(['GMRES iter RCIP = ',num2str(iter)])
%
% *** Construction of fine grid ***
  [sinterfin,sinterdifffin,npanfin]=panelinit(nse+iref,npan0);
  [zfin,~,nzfin,wzpfin,kappafin,awzpfin,zinterfin,npfin,sfin]= ...
      zinit(theta,sinterfin,sinterdifffin,T,W,npanfin);
%
% *** Reconstruction of layer density on the fine grid ***
  rhot=rhotilde(starind);
  pts2=Ng*(nse+2);
  rhofinloc=zeros(4*pts2,1);
  PbcL=Poper(IP);
  Pbc2L=blkdiag(PbcL,PbcL);
  for k=nse:-1:1
    Kcirc=Kstor(:,:,k);
    Kcirc(starL,starL)=zeros(8*Ng);
    MAT=eye(12*Ng)+Kcirc;
    tmp=Pbc2L*rhot;
    rhot=tmp-Kcirc*SchurBanaLs(MAT,Rstor(:,:,k),starL,circL,tmp,nse-k+1);
    myindL=[pts2-(k+2)*Ng+(1:Ng) pts2+(k+1)*Ng+(1:Ng)];
    rhofinloc([myindL 2*pts2+myindL])=rhot(circL);
    rhot=rhot(starL);
  end
  myindL=nse*Ng+1:(4+nse)*Ng;
  rhofinloc([myindL 2*pts2+myindL])=Rstor(:,:,1)*rhot;
%  clear Kstor Rstor
  rhofin=zeros(2*npfin,1);
  rhofin(1:pts2)=rhofinloc(pts2+1:2*pts2);
  rhofin(pts2+(1:np-4*Ng))=rhotilde(2*Ng+1:np-2*Ng);
  rhofin(npfin-pts2+1:npfin)=rhofinloc(1:pts2);
  rhofin(npfin+1:npfin+pts2)=rhofinloc(3*pts2+1:4*pts2);
  rhofin(npfin+pts2+(1:np-4*Ng))=rhotilde(np+2*Ng+1:2*np-2*Ng);
  rhofin(2*npfin-pts2+1:2*npfin)=rhofinloc(2*pts2+1:3*pts2);
  rho1fin=rhofin(1:npfin);
  rho2fin=rhofin(npfin+1:2*npfin);
%
% *** Zoom of layer densities close to corner ***
  if asymp==1
    figure(1)
% *** note that points 1:32 must be omitted. They are not rho, but rhohat
    myind=2*Ng+1:npfin/2;
    loglog(sfin(myind),abs(rho1fin(myind)),'b.', ...
	   sfin(myind),abs(rho2fin(myind)),'r.')
    title('Reconstruction of $\rho$ to the right of the corner', ...
	  'Interpreter','LaTeX')
    grid on
    xlabel('distance to corner','Interpreter','LaTeX')
    hl=legend('$|\rho_1|$','$|\rho_2|$','Location','NorthEast');
    set(hl,'Interpreter','LaTeX','FontSize',12)
    drawnow
  end
%
  rhohat=RG*rhotilde;
  rho1=rhohat(1:np);
  rho2=rhohat(np+1:2*np);
%
  [uref,uzref]=testfun(zg); % reference solution
  iout=Ineval(zg,theta,ngrtot);

  t1=cputime;
  [u,uz]=fieldcomp(rho1,rho2,z,nz,wzp,awzp,kappa,zinter,rho1fin,rho2fin, ...
      zfin,nzfin,wzpfin,awzpfin,kappafin,zinterfin,zg,dlim1,ngrtot,iout, ...
	  npan,npanfin,Ng,sfin,theta);
  t2=cputime-t1;
  disp(['The CPU time on evaluation is ',num2str(t2), ' seconds'])

%  fieldplot(u,uz,z,xg,yg,xylim,ngr)    
  errplot((1-iout).*u,(1-iout).*uref,(1-iout).*uz,(1-iout).*uzref, ...
	  xg,yg,xylim,ngr)
%
% find the maximum distance from points with large errors to the boundary
  ind=find(abs((1-iout).*(u-uref))>1e-13);
  ntmp=length(ind);
  disp(['# of points with potential error greater than 1e-13 = ',num2str(ntmp)])
  dmax=0;
  for i=ind'
    [d,~]=min(abs(zg(i)-zfin));
    if d>dmax
      dmax=d;
    end
  end
  if dmax>0
    figure(6)
    plot(real(z),imag(z),'b',real(zg(ind)),imag(zg(ind)),'r*', ...
	 real(zinter),imag(zinter),'g.')
    disp(['the maximum distance from these points to the boundary is ',num2str(dmax)])
  end

% find the maximum distance from points with large gradient errors to the boundary
  ind=find(abs((1-iout).*(uz-uzref))>1e-11);
  ntmp=length(ind);
  disp(['# of points with gradient error greater than 1e-11 = ',num2str(ntmp)])
  dmax=0;
  for i=ind'
    d=min(abs(zg(i)-zfin));
    if d>dmax
      dmax=d;
    end
  end
  if dmax >0
    figure(7)
    plot(real(z),imag(z),'b',real(zg(ind)),imag(zg(ind)),'r*', ...
	 real(zinter),imag(zinter),'g.')
    disp(['the maximum distance from these points to the boundary is ',num2str(dmax)])
  end 

  function [u,uz]=fieldcomp(rho1,rho2,z,nz,wzp,awzp,kappa,zinter,rho1fin, ...
	rho2fin,zfin,nzfin,wzpfin,awzpfin,kappafin,zinterfin,zg,dlim1, ...
	ngrtot,iout,npan,npanfin,Ng,sfin,theta)
  dlim4=0.04;
% *** Initial screening ***
  disp('Screening starts')
% *** 1 => regular evaluation
% *** 2 => some close evaluation
% *** 3 => some close corner-evaluation
  panlen=zeros(npan,1);
  for k=1:npan
    panlen(k)=sum(awzp((k-1)*Ng+1:k*Ng));
  end
  panlenfin=zeros(npanfin,1);
  for k=1:npanfin
    panlenfin(k)=sum(awzpfin((k-1)*Ng+1:k*Ng));
  end
  kvec=find(iout==0)';          % indicies of points inside
  mylist=zeros(ngrtot,1); 
  mylist(kvec)=1;
  t1=cputime;
  for k=kvec
    for kk=[1 npan]
      myind=(kk-1)*Ng+1:kk*Ng;
      d=min(abs(z(myind)-zg(k)))/panlen(kk);
      if d<2.0
	    mylist(k)=mylist(k)+2;  % close corner-evaluation needed
	    break
      end
    end     
    for kk=1:npan
      myind=(kk-1)*Ng+1:kk*Ng;
      d=min(abs(z(myind)-zg(k)))/panlen(kk);
      if d<dlim1 && mylist(k)==1
	    mylist(k)=mylist(k)+1;  % special evaluation needed
	    break
      end
    end

    for kk=1:npanfin
      myind=(kk-1)*Ng+1:kk*Ng;
      d=min(abs(zfin(myind)-zg(k)))/panlenfin(kk);
      if d<dlim4
	    mylist(k)=4;  % interpolation for very close evaluation needed    
      end
    end
  end

  mylist1=find(mylist==1)';
  mylist2=find(mylist==2)';
  mylist3=find(mylist==3)';
  mylist4=find(mylist==4)';
  evalpoints=length(mylist1)+length(mylist2)+length(mylist3) ...
	    +length(mylist4);
  disp(['Number of evaluation points = ',num2str(evalpoints)])

  [zclose,dclose,kclose]=findbdpoint(zg(mylist4),zfin,sfin, ...
      npanfin,panlenfin,theta,dlim4);
  
  u=zeros(ngrtot,1);
  uz=zeros(ngrtot,1);
  sorting_time=cputime-t1
%
  disp(['List1 starts, target points = ',num2str(length(mylist1))])
  for k=mylist1
    [u(k),uz(k)]=regueval(zg(k),z,nz,wzp,kappa,rho1,rho2);
  end
  disp(['List2 starts, target points = ',num2str(length(mylist2))])
  for k=mylist2
    [u(k),uz(k)]=speceval(zg(k),z,nz,wzp,kappa,zinter,rho1,rho2, ...
			  dlim1,panlen,npan,0,Ng);
  end
  disp(['List3 starts, target points = ',num2str(length(mylist3))])
  for k=mylist3
    [u(k),uz(k)]=speceval(zg(k),zfin,nzfin,wzpfin,kappafin, ...
        zinterfin,rho1fin,rho2fin,dlim1,panlenfin,npanfin,0,Ng);
  end
  disp(['List4 starts, target points = ',num2str(length(mylist4))])
  k4=1;
  for k=mylist4
    [u(k),uz(k)]=closeeval(zg(k),zclose(k4),dclose(k4),kclose(k4), ...
              zfin,nzfin,wzpfin,kappafin,zinterfin,...
              rho1fin,rho2fin,dlim1,panlenfin,npanfin,0,Ng);
    k4=k4+1;
  end
  
  u(iout==1)=mean(u(kvec)); % outside values set to something
%  
%
  function [zc,dc,kc]=findbdpoint(zt,z,s,npan,panlen,theta,dlim4)
% find the closed point zc, the shortest distance relative to the 
% panel length, and the corresponding panel to the given target
% points zt.
% Use Muller's method for root finding.
  nt=length(zt);
  zc=zeros(nt,1);
  dc=zeros(nt,1);
  kc=zeros(nt,1);
  nmax=100;
  myeps=1e-13;
  for k=1:nt
    for kk=1:npan
      myind=(kk-1)*16+(1:16);
      zs=z(myind);
      [d,id]=min(abs(zs-zt(k)));
      d=d/panlen(kk);
      if d<dlim4
        stmp=s(myind(id));
        kkmin=kk;
        break
      end
    end
    s1=stmp;
    s2=stmp+0.01;
    s3=stmp-0.01;
    [s0,~,~] = muller(nmax,myeps,theta,zt(k),s1,s2,s3);
    zc(k)=zfunc(s0,theta);
    dc(k)=abs(zt(k)-zc(k))/panlen(kkmin);
    kc(k)=kkmin;
  end    
%
  function [u,uz]=regueval(zt,z,nz,wzp,kappa,rho1,rho2)
  [K1,K2,K1z,K2z]=repkernels(zt,z,nz,wzp,kappa);
    
  u=(rho1.'*K1+rho2.'*K2)/pi;
  uz=(rho1.'*K1z+rho2.'*K2z)/pi;
   
  function [u,uz]=speceval(zt,z,nz,wzp,kappa,zinter,rho1,rho2,dlim1, ...
			   panlen,npan,io,Ng)
  [u,uz]=regueval(zt,z,nz,wzp,kappa,rho1,rho2);
  for k=1:npan
    myind=(k-1)*Ng+(1:Ng);
    zs=z(myind);
    d=min(abs(zs-zt))/panlen(k);
    if d<dlim1
      nzs=nz(myind);
      wzps=wzp(myind);
      kappas=kappa(myind);
      rho1s=rho1(myind);
      rho2s=rho2(myind);
      
      [CauC,HypC,SupC]=wCHSinit(zinter(k),zinter(k+1),zt,zs,wzps,io);
      [K1,K2,K1z,K2z]=repkernelsclose(CauC,HypC,SupC,zt,zs,nzs,kappas);
      u=u+(rho1s.'*K1+rho2s.'*K2)/pi;
      uz=uz+(rho1s.'*K1z+rho2s.'*K2z)/pi;
    end
  end

  function [u,uz]=closeeval(zt,zclose,dclose,kclose,z,nz,wzp, ...
			  kappa,zinter,rho1,rho2,dlim1,panlen,npan,io,Ng)
  n=8;
  x=(cos(pi*(n:-1:0)'/n)+1)/2; % Chebyshev nodes in ascending order
  if kclose<12 || kclose>npan-11
    x=x*1.0;
  else
    x=x*1.1;
  end
  zts=zclose+x*(zt-zclose)/dclose;  
  f=zeros(n+1,1);
  % value at the closest boundary point
  [~,f(1)]=testfun(zclose);
  for k=2:n+1
    [~,f(k)]=speceval(zts(k),z,nz,wzp,kappa,zinter,rho1,rho2,dlim1, ...
			   panlen,npan,io,Ng);
  end
  uz=chebbary(dclose,x,f);
  [u,~]=speceval(zt,z,nz,wzp,kappa,zinter,rho1,rho2,dlim1, ...
 			   panlen,npan,io,Ng);    
           
  function [K1,K2,K1z,K2z]=repkernels(ztarg,zs,nz,wzp,kappa)
  zm=zs-ztarg;
  rnu=real(conj(nz).*zm);
  f1=0.5*rnu.*nz;
  K1=real((f1./zm.^2+0.5./zm).*wzp/1i) ...
     -kappa.*real(conj(nz).*zm).*real(wzp./zm/1i)+0.5*kappa.*abs(wzp);
  K2=-0.5*rnu.*real(wzp./zm/1i)+0.25*abs(wzp);
  K1z=conj(2*f1.*wzp./zm.^3/1i+0.5*wzp./zm.^2/1i) ...
      -0.5*nz.*real(nz.*wzp./zm.^2/1i) ...
      +kappa.*(conj(-rnu.*wzp./zm.^2/1i)+nz.*real(wzp./zm/1i));
  K2z=conj(-0.5*rnu.*wzp./zm.^2/1i)+0.5*nz.*real(wzp./zm/1i);
    
  function [K1,K2,K1z,K2z]=repkernelsclose(CauC,HypC,SupC,ztarg,zs,nz,kappa)
  zm=zs-ztarg;
  f2=-0.5*real(conj(nz).*zm);
  f1=-f2.*nz;
  K1=real(f1.*HypC)+0.5*real(CauC)+2*f2.*kappa.*real(CauC);
  K2=f2.*real(CauC);
  K1z=conj(2*f1.*SupC+0.5*HypC)-0.5*nz.*real(nz.*HypC) ...
      +kappa.*(conj(2*f2.*HypC)+nz.*real(CauC));
  K2z=conj(f2.*HypC)+0.5*nz.*real(CauC);
 
  function fieldplot(u,uz,z,xg,yg,xylim,ngr)
  F1=zeros(ngr);
  for k=1:ngr
    F1(1:ngr,k)=u((k-1)*ngr+(1:ngr));
  end
  figure(2)
  imagesc(xg,yg,F1);      
  colormap(jet)
  axis xy
  colorbar
  hold on
  np=length(z)/2;
  xy=1.1*xylim;
  zext=[xy(1);0;z(1:np);1;xy(2);xy(2)+1i*xy(3);xy(1)+1i*xy(3);xy(1)];
  fill(real(zext),imag(zext),'w','EdgeColor','w')  
  zext=[xy(2);1;z(np+1:2*np);0;xy(1);xy(1)+1i*xy(4);xy(2)+1i*xy(4);xy(2)];
  fill(real(zext),imag(zext),'w','EdgeColor','w')  
  title('Field $u({\bf x})$','Interpreter','LaTeX')
  xlabel('$x_1$','Interpreter','LaTeX')
  ylabel('$x_2$','Interpreter','LaTeX')
  axis(xylim)
  axis equal
%
  F1=zeros(ngr);
  for k=1:ngr
    F1(1:ngr,k)=abs(uz((k-1)*ngr+(1:ngr)));
  end
  figure(3)
  imagesc(xg,yg,F1);      
  colormap(jet)
  axis xy
  colorbar
  hold on
  np=length(z)/2;
  xy=1.1*xylim;
  zext=[xy(1);0;z(1:np);1;xy(2);xy(2)+1i*xy(3);xy(1)+1i*xy(3);xy(1)];
  fill(real(zext),imag(zext),'w','EdgeColor','w')  
  zext=[xy(2);1;z(np+1:2*np);0;xy(1);xy(1)+1i*xy(4);xy(2)+1i*xy(4);xy(2)];
  fill(real(zext),imag(zext),'w','EdgeColor','w')  
  title('Field $|\nabla u({\bf x})|$','Interpreter','LaTeX')
  xlabel('$x_1$','Interpreter','LaTeX')
  ylabel('$x_2$','Interpreter','LaTeX')
  axis(xylim)
  axis equal
  drawnow

  function errplot(u,uref,uz,uzref,xg,yg,xylim,ngr)
  F1=zeros(ngr);
  for k=1:ngr
    F1(1:ngr,k)=abs(uref((k-1)*ngr+(1:ngr))-u((k-1)*ngr+(1:ngr)));
  end
  F1=log10(F1);
  maxa_err=max(max(F1));
  disp(['max field error = ',num2str(max(max(F1)))])
  set(groot,'DefaultAxesFontSize',15)
  figure(4)
  imagesc(xg,yg,F1,[log10(eps) maxa_err]);      
  colormap(flipud(pink(256)))
  axis xy
  colorbar
%  title('${\rm Log}_{10}$ of absolute error in $u({\bf x})$', ...
%	'Interpreter','LaTeX')
  xlabel('$x_1$','Interpreter','LaTeX')
  ylabel('$x_2$','Interpreter','LaTeX')
  axis(xylim)
  axis equal
%
  F1=zeros(ngr);
  for k=1:ngr
    F1(1:ngr,k)=abs(uzref((k-1)*ngr+(1:ngr))-uz((k-1)*ngr+(1:ngr)));
  end
  F1=log10(F1);
  maxa_err=max(max(F1));
  disp(['max gradient field error = ',num2str(max(max(F1)))])
  figure(5)
  imagesc(xg,yg,F1,[log10(eps) maxa_err]);      
  colormap(flipud(pink(256)))
  axis xy
  colorbar
%  title('${\rm Log}_{10}$ of absolute error in $|\nabla u({\bf x})|$', ...
%	'Interpreter','LaTeX')
  xlabel('$x_1$','Interpreter','LaTeX')
  ylabel('$x_2$','Interpreter','LaTeX')
  axis(xylim)
  axis equal
  
  function iout=Ineval(zg,theta,np)
  iout=zeros(np,1);            % point assumed to be inside
  thetag=abs(angle(zg));
  iout(thetag>theta/2)=1;     % point is outside
  sg=0.5+thetag/theta;
  iout(abs(zg)>sin(pi*sg))=1; % point is outside

  function [xg,yg,zg,ngrtot]=targinit(xylim,ngr)
  xg=linspace(xylim(1),xylim(2),ngr);
  yg=linspace(xylim(3),xylim(4),ngr);
  ngrtot=ngr^2;
  zg=zeros(ngrtot,1);
  for k=1:ngr
    zg((k-1)*ngr+(1:ngr))=xg(k)+1i*yg;
  end
  
  function [R,Rstor,Kstor]=Rcomp(R0,theta,T,W,Pbc,PWbc,nsub,prema, ...
				 nse,npan0,iref,starL,circL)
  Ng=length(T);
  Rstor=zeros(8*Ng,8*Ng,nse);
  Kstor=zeros(12*Ng,12*Ng,nse);
  myind1=[1:3*Ng 6*Ng+1:9*Ng];
  myind2=[3*Ng+1:6*Ng 9*Ng+1:12*Ng];
  denom=2^(nsub-1)*npan0*2^iref;
  D21=[ones(4*Ng) ones(4*Ng)/denom;ones(4*Ng)*denom ones(4*Ng)];
  R=D21.*R0;
% *** Recursion  
  for level=1:nsub
    [z,nz,wzp,kappa,awzp]=zlocinit(theta,T,W,nsub,level,npan0,iref);  
    MK=MKmatinit(z,nz,wzp,kappa,awzp,6*Ng);
    MK(myind1,myind1)=0;
    MK(myind2,myind2)=0;
    MAT=eye(12*Ng)+MK;
    if level >= prema
      Kstor(:,:,level-prema+1)=MK;
      Rstor(:,:,level-prema+1)=R;
    end
    R=SchurBana(Pbc,PWbc,MAT,R,starL,circL,nsub-level+1);
  end

  function R=initializeR0(T,W,theta,Pbc,PWbc,starL,circL)
  Ng=length(T);
  myind1=[1:3*Ng 6*Ng+1:9*Ng];
  myind2=[3*Ng+1:6*Ng 9*Ng+1:12*Ng];
  R=eye(8*Ng);
  myerr=1;
  iter=0;
  D21=[ones(4*Ng) 0.5*ones(4*Ng);2*ones(4*Ng) ones(4*Ng)];
  [z,nz,wzp,kappa,awzp]=zlocinit0(theta,T,W);  
  MK=MKmatinit(z,nz,wzp,kappa,awzp,6*Ng);
  MK(myind1,myind1)=0;
  MK(myind2,myind2)=0;
  MAT=eye(12*Ng)+MK;
  while myerr>2*eps && iter<1000
    Rold=R;
    R=D21.*SchurBana(Pbc,PWbc,MAT,R,starL,circL,1);
    myerr=norm(R-Rold,'fro')/norm(R,'fro');
    iter=iter+1;
  end
  disp(['Fixed-point iter  = ',num2str(iter)])
  disp(['Fixed-point error = ',num2str(myerr)])
  
  function A=SchurBana(P,PW,K,A,starL,circL,idepth)
  Ng=length(P)/8;
  starS=[Ng+1:3*Ng 5*Ng+1:7*Ng];
  circS=[1:Ng 3*Ng+1:5*Ng 7*Ng+1:8*Ng];
  VA=K(circL,starL)*A;
  PTA=PW'*A;
  PTAU=PTA*K(starL,circL);
  DVAUI=myinvnew(K(circL,circL)-VA*K(starL,circL),idepth);
  DVAUIVAP=DVAUI*(VA*P);
  A(starS,starS)=PTA*P+PTAU*DVAUIVAP;
  A(circS,circS)=DVAUI;
  A(circS,starS)=-DVAUIVAP;
  A(starS,circS)=-PTAU*DVAUI;

  function b=SchurBanaLs(K,A,starL,circL,rhs,idepth)
  VA=K(circL,starL)*A;
  AU=A*K(starL,circL);
  DVAU=K(circL,circL)-VA*K(starL,circL);
  b=zeros(length(rhs),1);
  DVAUIrhs=mylinsnew(DVAU,rhs(circL),idepth);
  DVAUIVArhs=mylinsnew(DVAU,VA*rhs(starL),idepth);
  b(starL)=A*rhs(starL)+AU*DVAUIVArhs-AU*DVAUIrhs;
  b(circL)=DVAUIrhs-DVAUIVArhs;
  
  function Mi=myinvnew(M,idepth)
% *** simple diagonal scaling ***
  np=length(M)/2;
  d1=[ones(np,1);ones(np,1)/2^idepth];
  d2=[ones(np,1);ones(np,1)*2^idepth];
%  koll1=cond(M)
%  koll2=cond(d1*d2'.*M)
%  pause
  Mi=d2(:,ones(1,2*np)).*((d1*d2'.*M)\diag(d1));
  
  function x=mylinsnew(M,b,idepth)
% *** simple diagonal scaling ***
  np=length(M)/2;
  d1=[ones(np,1);ones(np,1)/2^idepth];
  d2=[ones(np,1);ones(np,1)*2^idepth];
  x=d2.*((d1*d2'.*M)\(d1.*b));  
  
  function [x,it]=myGMRESR(A,R,b,n,m,tol)
% *** GMRES with low-threshold stagnation control ***
  V=zeros(n,m+1);
  H=zeros(m);
  cs=zeros(m,1);
  sn=zeros(m,1);
  bnrm2=norm(b);
  V(:,1)=b/bnrm2;
  s=bnrm2*eye(m+1,1);
  for it=1:m                                  
    it1=it+1;                                   
    w=A*(R*V(:,it));
    for k=1:it
      H(k,it)=V(:,k)'*w;
      w=w-H(k,it)*V(:,k);
    end
    H(it,it)=H(it,it)+1;
    wnrm2=norm(w);
    V(:,it1)=w/wnrm2;
    for k=1:it-1                                
      temp     = cs(k)*H(k,it)+sn(k)*H(k+1,it);
      H(k+1,it)=-sn(k)*H(k,it)+cs(k)*H(k+1,it);
      H(k,it)  = temp;
    end
    [cs(it),sn(it)]=rotmat(H(it,it),wnrm2);     
    H(it,it)= cs(it)*H(it,it)+sn(it)*wnrm2;
    s(it1) =-sn(it)*s(it);                      
    s(it)  = cs(it)*s(it);                         
    myerr=abs(s(it1))/bnrm2;
    if (myerr<=tol)||(it==m)                     
      disp(['predicted residual = ' num2str(myerr)])
      y=triu(H(1:it,1:it))\s(1:it);             
      x=fliplr(V(:,1:it))*flipud(y);
      trueres=norm(x+A*(R*x)-b)/bnrm2;
      disp(['true residual      = ',num2str(trueres)])                  
      break
    end
  end
  
  function [c,s]=rotmat(a,b)
  if  b==0
    c=1;
    s=0;
  elseif abs(b)>abs(a)
    temp=a/b;
    s=1/sqrt(1+temp^2);
    c=temp*s;
  else
    temp=b/a;
    c=1/sqrt(1+temp^2);
    s=temp*c;
  end

  function MK=MKmatinit(z,nz,wzp,kappa,awzp,N)
% *** 4 kernels for 2x2 integral equation system ***
  K11=zeros(N); K12=zeros(N); K21=zeros(N); K22=zeros(N);
  for m=1:N
    zm=z(m)-z;
    zm0=conj(zm)./zm;
    K11(:,m)= 0.25*imag((3+nz(m)^2*zm0)*wzp(m)./zm) ...
	     -0.5*kappa(m)*imag(nz(m)*wzp(m)*zm0);	      
    K12(:,m)=-0.25*imag(nz(m)*wzp(m)*zm0);
    K21(:,m)= 0.25*imag((3*nz-nz(m)^2*(conj(nz)-2*nz.*zm0))*wzp(m)./zm.^2) ...
	     +0.5*kappa(m)*imag(nz(m)*(conj(nz)-nz.*zm0)*wzp(m)./zm);
    K22(:,m)= 0.25*imag(nz(m)*(conj(nz)-nz.*zm0)*wzp(m)./zm);
  end
% *** correct diagonal terms ***
  K11(1:N+1:end)=0.5*kappa.*awzp;
  K12(1:N+1:end)=0.25*awzp;
  K21(1:N+1:end)=0.25*kappa.^2.*awzp;
  K22(1:N+1:end)=0.5*kappa.*awzp;
  MK=2/pi*[K11 K12; K21 K22]; % assemble the matrix blocks
  
  function [z,zp,nz,wzp,kappa,awzp,zinter,np,s]=zinit(theta,sinter, ...
						      sinterdiff,T,W,npan)
  Ng=length(T);
  np=Ng*npan;
  s=zeros(np,1);
  w=zeros(np,1);
  for k=1:npan
    myind=(k-1)*Ng+1:k*Ng;
    sdif=sinterdiff(k)/2;
    s(myind)=(sinter(k)+sinter(k+1))/2+sdif*T;
    w(myind)=W*sdif;
  end
  z=zfunc(s,theta) ;
  zp=zpfunc(s,theta);
  zpp=zppfunc(s,theta);
  zinter=zfunc(sinter,theta);
% *** some extra presicion gained from symmetry ***
  z(np/2+1:np)=conj(flipud(z(1:np/2)));
  zp(np/2+1:np)=-conj(flipud(zp(1:np/2)));
  zpp(np/2+1:np)=conj(flipud(zpp(1:np/2)));
  itmp=round(npan/2+1);
  zinter(itmp:npan+1)=conj(flipud(zinter(1:npan-itmp+2)));
% *************************************************
  nz=-1i*zp./abs(zp);
  wzp=w.*zp;
  awzp=w.*abs(zp);
  kappa = imag(zpp./zp)./abs(zp);

  function [z,nz,wzp,kappa,awzp]=zlocinit0(theta,T,W)
  Ng=length(T);
  s=[T/4+0.25;T/4+0.75;T/2+1.5];
  w=[W/4;W/4;W/2];
  w=[flipud(w);w];
  z=pi*s.*exp(-1i*theta*0.5);
  z=[conj(flipud(z));z];
  zp=pi*exp(-1i*theta*0.5)*ones(3*Ng,1);
  zp=[-conj(flipud(zp));zp];
  nz=-1i*zp./abs(zp);
  wzp=w.*zp;
  awzp=w.*abs(zp);
  kappa=zeros(6*Ng,1);
  
  function [z,nz,wzp,kappa,awzp]=zlocinit(theta,T,W,nsub,level,npan0,iref)
  denom=2^(nsub-level)*npan0*2^iref;
  s=[T/4+0.25;T/4+0.75;T/2+1.5]/denom;
  w=[W/4;W/4;W/2]/denom;
  w=[flipud(w);w];
  z=zfunc(s,theta);  
  z=[conj(flipud(z));z];
  zp=zpfunc(s,theta);
  zp=[-conj(flipud(zp));zp];
  zpp=zppfunc(s,theta);
  zpp=[conj(flipud(zpp));zpp];
  nz=-1i*zp./abs(zp);
  wzp=w.*zp;
  awzp=w.*abs(zp);
  kappa = imag(zpp./zp)./abs(zp);
  
  function zout=zfunc(s,theta)
  zout=sin(pi*s).*exp(1i*theta*(s-0.5));

  function zpout=zpfunc(s,theta)
  zpout=(pi*cos(pi*s)+1i*theta*sin(pi*s)).*exp(1i*theta*(s-0.5));
  
  function zppout=zppfunc(s,theta)
  zppout=(2i*pi*theta*cos(pi*s)-(theta^2+pi^2)*sin(pi*s)).* ...
	 exp(1i*theta*(s-0.5));  

  function [sinter,sinterdiff,npanfin]=panelinit(nsub,npan)
  npanfin=npan+2*nsub;
  sinter=zeros(npanfin+1,1);
  sinter(1:npan+1)=linspace(0,1,npan+1);
  sinterdiff=ones(npan+2*nsub,1)/npan;
  for k=1:nsub
    sinter(3:end)=sinter(2:end-1);
    sinter(2)=(sinter(1)+sinter(2))/2;   
    sinterdiff(3:end)=sinterdiff(2:end-1);
    sinterdiff(2)=sinterdiff(1)/2;
    sinterdiff(1)=sinterdiff(1)/2;   
  end
  sinter(end-nsub:end)=1-flipud(sinter(1:nsub+1));
  sinterdiff(end-nsub-1:end)=flipud(sinterdiff(1:nsub+2));
  
  function [wcmpC,wcmpH,wcmpS]=wCHSinit(a,b,ztg,zsc,wzp,io)
% *** ztgtr is target vector in transformed plane ***
  Ng=length(zsc);
  cc=(b-a)/2;
  ztgtr=(ztg-(b+a)/2)/cc;
  zsctr=(zsc-(b+a)/2)/cc;
  A=ones(Ng);
  for k=2:Ng
    A(:,k)=zsctr.*A(:,k-1);
  end
  p=zeros(Ng,1);
  r=zeros(Ng,1);
  s=zeros(Ng,1);
  c1=(1-(-1).^(1:Ng))./(1:Ng);
  c2=-(1/(1-ztgtr)+(-1).^(0:Ng-1)/(1+ztgtr));
  c3=-0.5*(1/(1-ztgtr)^2+(-1).^(1:Ng)/(1+ztgtr)^2);
  upp=log(1-ztgtr);
  loo=log(-1-ztgtr);
  if io==1 % point is outside
    if imag(ztgtr)>0 && abs(real(ztgtr))<1
      loo=loo+2i*pi;
    end
  else     % point is inside
    if imag(ztgtr)<0 && abs(real(ztgtr))<1
      loo=loo-2i*pi;    
    end
  end
  p(1)=upp-loo;
  r(1)=c2(1);
  s(1)=c3(1);
  for k=1:Ng-1
    p(k+1)=ztgtr*p(k)+c1(k);
    r(k+1)=k*p(k)+c2(k+1);  
    s(k+1)=0.5*k*r(k)+c3(k+1);
  end
  wcmpC=A.'\p/1i-wzp./(zsc-ztg)/1i;
  wcmpH=A.'\r/1i/cc-wzp./(zsc-ztg).^2/1i;
  wcmpS=A.'\s/1i/cc.^2-wzp./(zsc-ztg).^3/1i;
  
  function [IP,IPW]=IPinit(T,W)
  Ng=length(T);
  A=ones(Ng);
  AA=ones(2*Ng,Ng);
  T2=[T-1;T+1]/2;
  W2=[W;W]/2;
  for k=2:Ng
    A(:,k)=A(:,k-1).*T;
    AA(:,k)=AA(:,k-1).*T2;   
  end
  IP=AA/A;
  IPW=IP.*(W2*(1./W)');

  function Pbc=Pbcinit(IP)
  Pbc=blkdiag(IP,rot90(IP,2));

  function Pbc=Poper(IP)
  Ng=length(IP)/2;
  Pbc=blkdiag(eye(Ng),IP,rot90(IP,2),eye(Ng));

  function T=Tinit16
% *** 16-point Gauss-Legendre nodes ***  
  T=zeros(16,1);
  T( 1)=-0.989400934991649932596154173450332627;
  T( 2)=-0.944575023073232576077988415534608345;
  T( 3)=-0.865631202387831743880467897712393132;
  T( 4)=-0.755404408355003033895101194847442268;
  T( 5)=-0.617876244402643748446671764048791019;
  T( 6)=-0.458016777657227386342419442983577574;
  T( 7)=-0.281603550779258913230460501460496106;
  T( 8)=-0.095012509837637440185319335424958063;
  T( 9)= 0.095012509837637440185319335424958063;
  T(10)= 0.281603550779258913230460501460496106;
  T(11)= 0.458016777657227386342419442983577574;
  T(12)= 0.617876244402643748446671764048791019;
  T(13)= 0.755404408355003033895101194847442268;
  T(14)= 0.865631202387831743880467897712393132;
  T(15)= 0.944575023073232576077988415534608345;
  T(16)= 0.989400934991649932596154173450332627;

  function W=Winit16
% *** 16-point Gauss-Legendre weights ***  
  W=zeros(16,1); 
  W( 1)= 0.027152459411754094851780572456018104;
  W( 2)= 0.062253523938647892862843836994377694;
  W( 3)= 0.095158511682492784809925107602246226;
  W( 4)= 0.124628971255533872052476282192016420;
  W( 5)= 0.149595988816576732081501730547478549;
  W( 6)= 0.169156519395002538189312079030359962;
  W( 7)= 0.182603415044923588866763667969219939;
  W( 8)= 0.189450610455068496285396723208283105;
  W( 9)= 0.189450610455068496285396723208283105;
  W(10)= 0.182603415044923588866763667969219939;
  W(11)= 0.169156519395002538189312079030359962;
  W(12)= 0.149595988816576732081501730547478549;
  W(13)= 0.124628971255533872052476282192016420;
  W(14)= 0.095158511682492784809925107602246226;
  W(15)= 0.062253523938647892862843836994377694;
  W(16)= 0.027152459411754094851780572456018104;
  
  function [f,fz,zs]=testfun(ztarg)
% *** the fundamental solution of the modified biharmonic equation in 2D ***
  q=[0.171;0.720; 0.918; 0.825; 0.078];
  x=[1.280;0.413;-0.235;-0.145; 1.050];
  y=[0.177;0.795; 0.315;-0.474;-0.580];
  zs=x+1i*y;
  ns=5;
  f = zeros(size(ztarg));
  if nargout >1 
    fz = zeros(size(ztarg)); 
  end
  for i=1:ns
    z=ztarg-zs(i);
    r=abs(z);
    f = f+q(i)*r.^2.*log(r);
    if nargout > 1
      fz =fz+q(i)*z.*(2*log(r)+1);
    end
  end

  function [z,niter,err] = muller(nmax,eps,theta,zt,z1,z2,z3)
% Muller's method for root finding
% taylored for finding a point on the curve whose normal vector
% passes through the target point zt.
%
% nmax - maximum number of iterations
% eps - desired precision
% theta - the angle parameter for the curve
% zt - the target point
% z1,z2,z3 - three initial guesses for the root
%
% initialize variables
  niter = 0;                      % counts iteration steps
  err=100*eps;                    % percent error
  if nargin == 4
    x = rand(1,3);   % 3 initial guesses
    z1 = x(1); z2 = x(2); z3 = x(3);
  end

  w1 = mullerfunc(z1,theta,zt);
  w2 = mullerfunc(z2,theta,zt);
  w3 = mullerfunc(z3,theta,zt);

% iterate until max iterations or tolerance is met
  while niter < nmax && err>eps
    niter = niter + 1  ;        % update iteration step
    
    h1=z2-z1;
    h2=z3-z2;
    lambda=h2/h1;
    g=w1*lambda*lambda-w2*(1+lambda)*(1+lambda)+w3*(1+2*lambda);
    det=g*g-4*w3*(1+lambda)*lambda*(w1*lambda-w2*(1+lambda)+w3);
    
    h1=g+sqrt(det);
    h2=g-sqrt(det);
    %    
    if (abs(h1)>abs(h2))
        lambda=-2*w3*(1+lambda)/h1;
    else
        lambda=-2*w3*(1+lambda)/h2;
    end
    
    z1=z2;
    w1 = w2;
    z2=z3;
    w2 = w3;
    z3 = z2+lambda*(z2-z1);
    w3 = mullerfunc(z3,theta,zt);

    err=abs((z3-z2)/z3);    
  end
  z=z3;
%
%
%
  function f = mullerfunc(x,theta,zt)
% return the value of the Legendre expansion of order n at x,
% the expansion coefficients are given by coef.
  z=zfunc(x,theta);
  zp=zpfunc(x,theta);
  zn=-1i*zp./abs(zp);
  f=imag((zt-z)*conj(zn));
%
%
%
  function ff=chebbary(xx,x,f)
% barycentric interpolation using Chebyshev nodes of the second kind,
% the nodes x could be scaled and shifted. 
% Nodes are in the ascending order!
%
% From Jean-Paul Berrut and Lloyd N. Trefethen
% "Barycentric Lagrange Interpolation" SIAM Review, Vol. 46, No. 3, 
% pp. 501–517, 2004.
%
  n=length(x)-1;
% c contains Chebyshev interpolation weights. It
% remain the same if nodes are shifted and/or scaled
% due to cancellation of common factors.
%
% It can be computed in O(n) operations.
% Equispaced nodes is the only other case where 
% the weights can be computed in O(n) time.
  c=[1/2; ones(n-1,1); 1/2].*(-1).^((n:-1:0)'); 
  numer = zeros(size(xx));
  denom = zeros(size(xx));
  exact = zeros(size(xx));
  for j = 1:n+1
    xdiff = xx-x(j);
    temp = c(j)./xdiff;
    numer = numer + temp*f(j);
    denom = denom + temp;
    exact(xdiff==0) = 1;
  end
  ff = numer./denom;
  jj = find(exact); ff(jj) = f(exact(jj));


  
