  function modex3j
% *****************************************************************
% * Produces bottom row in Figure 8 in                            *
% * "On integral equation methods for the first Dirichlet problem *
% *  of the biharmonic and modified biharmonic equations"         *
% * SIAM J. Sci. Comput., vol. 40(4), pp. A2609-A2630, 2018.      *
% * https://doi.org/10.1137/17M1162238                            *
% *****************************************************************
  close all
  format short E
  format compact
% *** Comments on the parameters ************************************
%
% iref:  Coarse mesh refinement towards corner vertex. Recommended
%        value is iref=0. A larger iref, in combination with dlim1=1,
%        could lead to speedups in the post-processor.
%
% *** User specified quantities *************************************
  lambda=2.0e1     % parameter in the modified biharmonic equation
  theta=1.0*pi/2;  % corner opening angle
  t=theta*2/pi;    % tested case: t \in [0.4,1.7]
  if t<0.5
    npan0=42       % number of un-refined coarse panels
  elseif t<0.7
    npan0=36
  elseif t<1
    npan0=32
  else
    npan0=45
  end
  iref=0;        % Coarse mesh refinement towards corner vertex
  nsub=40;       % number of subdivisions
  dlim1=0.85;    % closeness parameter for 16-point quadrature
  ngr=200;       % field evaluation at ngr^2 points
  xylim=[-0.1 1.1 -0.5353 0.5353];  % computational domain
% **************************************************
  T=Tinit16;
  W=Winit16;
  [IP,IPW]=IPinit(T,W);
  Pbc=Pbcinit(IP);
  PWbc=Pbcinit(IPW);
  Pbc2=blkdiag(Pbc,Pbc);
  PWbc2=blkdiag(PWbc,PWbc);
  sidiloc=[1;0.5;0.5;0.5;0.5;1];
  LogCloc=LogCinit(sidiloc,T,W,6,96,0);
%
% *** Panels, points, and weights ***
  [sinter,sinterdiff,npan]=panelinit(iref,npan0);
  [z,zp,nz,wzp,kappa,awzp,zinter,np]=zinit(theta,sinter,sinterdiff,T,W,npan);
  ntot=2*np;
  [xg,yg,zg,ngrtot]=targinit(xylim,ngr);
%
% *** The whole interaction matrix M_coa^{\circ} is set up ***
  LogC=LogCinit(sinterdiff,T,W,npan,np,1);
  MS=MSmatinit(z,nz,kappa,awzp,lambda,np); %  smooth part of MAT
  ML=MLmatinit(LogC,z,zp,nz,kappa,awzp,sinterdiff,lambda,np); % log part of MAT
  starind=[np-31:np 1:32 2*np-31:2*np np+1:np+32];
% set all nonsmooth blocks of MS to zero
  myind1=[np-31:np 2*np-31:2*np ];
  myind2=[1:32 np+1:np+32];
  MS(myind1,myind2)=zeros(64);
  MS(myind2,myind1)=zeros(64);  
  ML(starind,starind)=zeros(128);
  Mcirc=MS+ML;  
%
% *** Construct the right hand side ***
  [g,gz]=testfun(z,lambda);
  gn=real(conj(nz).*gz);
  rhs=[g;gn];  
%  
% *** Recursion for the R matrix ***
  starL=[17:80 113:176];
  circL=[1:16 81:112 177:192];
  prema=min(nsub-ceil(log2(sum(awzp(1:16))/min(abs(zg)))),nsub);
  disp(['nsub = ',num2str(nsub)])
  %prema=min(nsub-ceil(log2(sum(awzp(1:16))/min(abs(zg)))),nsub);
  nse=nsub-prema+1;
  disp(['nse = ',num2str(nse)])
  [R,Rstor,Kstor]=Rcomp(theta,lambda,LogCloc,T,W,Pbc2,PWbc2,...
		      nsub,prema,nse,npan0,iref,starL,circL);
  RG=speye(ntot);
  RG(starind,starind)=R;  
%
% *** Solving main linear system ***
  [rhotilde,iter,~]=myGMRESR(Mcirc,RG,rhs,ntot,100,eps);
  disp(['GMRES iter RCIP = ',num2str(iter)]) 
%
% *** Construction of fine grid ***
  [sinterfin,sinterdifffin,npanfin]=panelinit(nse+iref,npan0);
  [zfin,~,nzfin,wzpfin,kappafin,awzpfin,zinterfin,npfin,sfin]= ...
      zinit(theta,sinterfin,sinterdifffin,T,W,npanfin);
%
% *** Reconstruction of layer density on the fine grid ***
  rhot=rhotilde(starind);
  pts2=16*(nse+2);
  rhofinloc=zeros(4*pts2,1);
  PbcL=Poper(IP);
  Pbc2L=blkdiag(PbcL,PbcL);
  for k=nse:-1:1
    Kcirc=Kstor(:,:,k);
    Kcirc(starL,starL)=zeros(128);
    MAT=eye(192)+Kcirc;
    tmp=Pbc2L*rhot;
    rhot=tmp-Kcirc*SchurBanaLs(MAT,Rstor(:,:,k),starL,circL,tmp);
    myindL=[pts2-(k+2)*16+(1:16) pts2+(k+1)*16+(1:16)];
    rhofinloc([myindL 2*pts2+myindL])=rhot(circL);
    rhot=rhot(starL);
  end
  myindL=nse*16+1:(4+nse)*16;
  rhofinloc([myindL 2*pts2+myindL])=Rstor(:,:,1)*rhot;
  rhofin=zeros(2*npfin,1);
  rhofin(1:pts2)=rhofinloc(pts2+1:2*pts2);
  rhofin(pts2+(1:np-64))=rhotilde(33:np-32);
  rhofin(npfin-pts2+1:npfin)=rhofinloc(1:pts2);
  rhofin(npfin+1:npfin+pts2)=rhofinloc(3*pts2+1:4*pts2);
  rhofin(npfin+pts2+(1:np-64))=rhotilde(np+33:2*np-32);
  rhofin(2*npfin-pts2+1:2*npfin)=rhofinloc(2*pts2+1:3*pts2);
  rho1fin=rhofin(1:npfin);
  rho2fin=rhofin(npfin+1:2*npfin);
%
  rhohat=RG*rhotilde;
  rho1=rhohat(1:np);
  rho2=rhohat(np+1:2*np);
%
  [uref,uzref]=testfun(zg,lambda); % reference solution
  iout=Ineval(zg,theta,ngrtot);
  [u,uz]=fieldcomp(lambda,rho1,rho2,z,nz,wzp,kappa,awzp,zinter,rho1fin, ...
   rho2fin,zfin,nzfin,wzpfin,kappafin,awzpfin,zinterfin,zg,dlim1,ngrtot,...
   iout,npan,npanfin,sfin,theta);
  %fieldplot(u,uz,z,xg,yg,xylim,ngr)
  errplot((1-iout).*u,(1-iout).*uref,(1-iout).*uz,(1-iout).*uzref, ...
	  xg,yg,xylim,ngr)
%
% find the maximum distance from points with large errors to the boundary
  ind=find(abs((1-iout).*(u-uref))>1e-13);
  ntmp=length(ind);
  disp(['# of points with potential error greater than 1e-13 = ',num2str(ntmp)])
  dmax=0;
  for i=ind'
    d=min(abs(zg(i)-zfin));
    if d>dmax
      dmax=d;
    end
  end
  if dmax>0
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
  if dmax>0
    disp(['the maximum distance from these points to the boundary is ',num2str(dmax)])
  end 

  function [u,uz]=fieldcomp(lambda,rho1,rho2,z,nz,wzp,kappa,awzp,zinter, ...
	rho1fin,rho2fin,zfin,nzfin,wzpfin,kappafin,awzpfin,zinterfin,zg, ...
	dlim1,ngrtot,iout,npan,npanfin,sfin,theta)
  dlim4=0.04;
% *** Initial screening ***
  disp('Screening starts')
% *** 1 => regular evaluation
% *** 2 => some close evaluation
% *** 3 => some close corner-evaluation
  panlen=zeros(npan,1);
  for k=1:npan
    panlen(k)=sum(awzp((k-1)*16+1:k*16));
  end
  panlenfin=zeros(npanfin,1);
  for k=1:npanfin
    panlenfin(k)=sum(awzpfin((k-1)*16+1:k*16));
  end
  kvec=find(iout==0)';          % indicies of points inside
  mylist=zeros(ngrtot,1); 
  mylist(kvec)=1;
  for k=kvec
    for kk=[1 npan]
      myind=(kk-1)*16+1:kk*16;
      d=min(abs(z(myind)-zg(k)))/panlen(kk);
      if d<2.8
    	mylist(k)=mylist(k)+2;  % close corner-evaluation needed
	break
      end
    end
    for kk=2:npan-1
      myind=(kk-1)*16+1:kk*16;
      d=min(abs(z(myind)-zg(k)))/panlen(kk);
      if d<dlim1 && mylist(k)==1
	mylist(k)=mylist(k)+1;  % special evaluation needed
	break
      end
    end
    for kk=1:npanfin
      myind=(kk-1)*16+(1:16);
      zs=zfin(myind);
      d=min(abs(zs-zg(k)))/panlenfin(kk);
      if d<dlim4
        mylist(k)=4; % interpolation used for close evaluation
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
%
  disp(['List1 starts, target points = ',num2str(length(mylist1))])
  for k=mylist1
    [u(k),uz(k)]=regueval(lambda,zg(k),z,nz,kappa,awzp,rho1,rho2);   
  end
  disp(['List2 starts, target points = ',num2str(length(mylist2))])
  for k=mylist2
    [u(k),uz(k)]=speceval(lambda,zg(k),z,nz,wzp,kappa,awzp,zinter,rho1,rho2,...
			  dlim1,panlen,npan,0);
  end
  disp(['List3 starts, target points = ',num2str(length(mylist3))])
  for k=mylist3
    [u(k),uz(k)]=speceval(lambda,zg(k),zfin,nzfin,wzpfin,kappafin,awzpfin, ...
			  zinterfin,rho1fin,rho2fin,dlim1,panlenfin,npanfin,0);
  end
  disp(['List4 starts, target points = ',num2str(length(mylist4))])
  k4=1;
  for k=mylist4
    [u(k),uz(k)]=closeeval(lambda,zg(k),zclose(k4),dclose(k4),kclose(k4), ...
        zfin,nzfin,wzpfin,kappafin,awzpfin,zinterfin, ...
        rho1fin,rho2fin,dlim1,panlenfin,npanfin,0);
    k4=k4+1;
  end 
  u(iout==1)=mean(u(kvec)); % outside values set to something

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

  function [u,uz]=regueval(lambda,zt,z,nz,kappa,awzp,rho1,rho2)
  [K1,K2,K1z,K2z]=repkernels(zt,z,nz,kappa,lambda);
  u=awzp'*(rho1.*K1+rho2.*K2);
  u=u*2/pi;

  uz=awzp'*(rho1.*K1z+rho2.*K2z);
  uz=uz*2/pi;

  function [u,uz]=speceval(lambda,zt,z,nz,wzp,kappa,awzp,zinter,rho1,rho2,...
                       dlim1,panlen,npan,io)
  [u,uz]=regueval(lambda,zt,z,nz,kappa,awzp,rho1,rho2);
  for k=1:npan
    myind=(k-1)*16+(1:16);
    zs=z(myind);
    d=min(abs(zs-zt))/panlen(k);    
    if d<dlim1
      nzs=nz(myind);
      wzps=wzp(myind);
      kappas=kappa(myind);
      awzps=awzp(myind);
      rho1s=rho1(myind);
      rho2s=rho2(myind);

      [LogC,CauC,HypC,SupC]=wLCHSinit(zinter(k),zinter(k+1),zt,zs,...
                                     nzs,wzps,awzps,io);
      [K1C,K2C,K1zC,K2zC]=repkernelsclose(LogC,CauC,HypC,SupC,zt,zs,...
					  nzs,kappas,awzps,lambda);
      u=u+(rho1s.'*K1C+rho2s.'*K2C)*2/pi;
      uz=uz+(rho1s.'*K1zC+rho2s.'*K2zC)*2/pi;
    end
  end

  function [u,uz]=closeeval(lambda,zt,zclose,dclose,kclose, ...
          z,nz,wzp,kappa,awzp,zinter,rho1,rho2,dlim1, ...
          panlen,npan,io)
  n=8;
  x=(cos(pi*(n:-1:0)'/n)+1)/2; % Chebyshev nodes in ascending order  
  if kclose<12 || kclose>npan-11
    x=x*1.0;
  else
    x=x*1.1;
  end
  zts=zclose+x*(zt-zclose)/dclose;
  f=zeros(n+1,1);
  [~,f(1)]=testfun(zclose,lambda); % value at the closest boundary point
  for k=2:n+1
    [~,f(k)]=speceval(lambda,zts(k),z,nz,wzp,kappa,awzp,zinter, ...
                      rho1,rho2,dlim1,panlen,npan,io);
  end
  uz=chebbary(dclose,x,f);
  [u,~]=speceval(lambda,zt,z,nz,wzp,kappa,awzp,zinter,rho1,rho2,dlim1, ...
 			   panlen,npan,io);

  function [K1,K2,K1z,K2z]=repkernels(ztarg,zs,nz,kappa,lambda)
% kernels for the representation of the solution and its derivatives
  zm = ztarg-zs;
  r1=abs(zm);
  r2=r1.*r1;
  r4=r2.*r2;
  rnu=real(nz.*conj(zm));
  rnu2=rnu.*rnu;
  rnu3=rnu.*rnu2;
  z=lambda*r1;
  z2=z.^2;
  z4=z2.*z2;
  [c0,c1,c2,c3,c4,c5]=cjseval(z);
  K1=-rnu./r2.*c1+rnu3./r4.*(c2-1+z2/8);
  K2=(-0.5+rnu2./r2).*(c0-0.5);
  K1=K1+2*kappa.*K2;
  K1z = (-nz.*c1+3*nz.*rnu2./r2.*(c2-1+z2/8) + zm.*rnu./r2.*(c3-z2/8) ...
	 -zm.*rnu3./r4.*(c4-4+z2/4-z4/48))./r2;    
  K2z = (2*rnu.*nz.*(c0-0.5) + 0.5*zm.*c5 -zm.*rnu2./r2.*(c2-1+z2/8))./r2;
  K1z=K1z+2*kappa.*K2z;

  function [K1C,K2C,K1zC,K2zC]=repkernelsclose(LogC,CauC,HypC,SupC,ztarg,zs,...
					       nz,kappa,awzp,lambda)
% corrections for near-field evaluation
  zm = ztarg-zs;
  r1=abs(zm);
  r2=r1.*r1;
  r4=r2.*r2;
  r6=r2.*r4;
  rnu=real(conj(nz).*zm);
  rnu2=rnu.*rnu;
  rnu3=rnu.*rnu2;
  z=lambda*r1;
  % Now find the logarithmic part of the kernels
  [c0lg,c1lg,c2lg,c3lg,c4lg,c5lg]=cjsLOGeval(z);
  c0lg=c0lg./r2;
  c1lg=c1lg./r2;
  c2lg=c2lg./r4;
  c3lg=c3lg./r4;
  c4lg=c4lg./r6;
  c5lg=c5lg./r2;
  f1a=-0.5*rnu.*nz;
  %K1=-rnu./r2.*c1+rnu3./r4.*c2;
  K1C=(-rnu.*c1lg+rnu3.*c2lg).*LogC.*awzp ...
      +real(f1a.*HypC)+(0.5-lambda^2/8*rnu2).*real(CauC);
  %K2=(-0.5+rnu2./r2).*c0;
  K2C=(-0.5*r2+rnu2).*c0lg.*LogC.*awzp + 0.5*rnu.*real(CauC);
  K1C=K1C+2*kappa.*K2C;
  %K1z = (-nz.*c1+3*nz.*rnu2./r2.*c2 + zm.*rnu./r2.*c3 -zm.*rnu3./r4.*c4)./r2;
  K1zlg = -nz.*c1lg+3*nz.*rnu2.*c2lg + zm.*rnu.*c3lg -zm.*rnu3.*c4lg;    
  %K2z = (2*rnu.*nz.*c0 + 0.5*zm.*c5 -zm.*rnu2./r2.*c2)./r2;
  K2zlg = 2*rnu.*nz.*c0lg + 0.5*zm.*c5lg -zm.*rnu2.*c2lg;
  K1zC=K1zlg.*LogC.*awzp...
      -lambda^2/8*(-2*zm+3*nz.*rnu+1/6*lambda^2*zm.*rnu2).*real(CauC)...
      +0.5*conj(HypC)-(0.5*nz+lambda^2/8*zm.*rnu).*real(nz.*HypC)...    
      -conj(rnu.*nz.*SupC);
  K2zC=K2zlg.*LogC.*awzp+0.5*(nz+0.25*lambda^2*zm.*rnu).*real(CauC)...
      +0.5*rnu.*conj(HypC);
  K1zC=K1zC+2*kappa.*K2zC;

  function fieldplot(u,uz,z,xg,yg,xylim,ngr)
  set(0,'DefaultAxesFontSize',12)
  F1=zeros(ngr);
  for k=1:ngr
    F1(1:ngr,k)=u((k-1)*ngr+(1:ngr));
  end
  fh = findobj('type','figure');
  nf = length(fh);
  figure(nf+1)
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
  fh = findobj('type','figure');
  nf = length(fh);
  figure(nf+1)
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
  set(0,'DefaultAxesFontSize',12)
  F1=zeros(ngr);
  for k=1:ngr
    F1(1:ngr,k)=abs(uref((k-1)*ngr+(1:ngr))-u((k-1)*ngr+(1:ngr)));
  end
  F1=log10(F1);
  maxa_err=max(max(F1));
  disp(['max field error = ',num2str(max(max(F1)))])
  fh =  findobj('type','figure');
  nf = length(fh);
  figure(nf+1)
  imagesc(xg,yg,F1,[log10(eps) maxa_err]);      
  colormap(flipud(pink(256)))
  axis xy
  colorbar
  title('${\rm Log}_{10}$ of absolute error in $u({\bf x})$', ...
	'Interpreter','LaTeX')
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
  fh =  findobj('type','figure');
  nf = length(fh);
  figure(nf+1)
  imagesc(xg,yg,F1,[log10(eps) maxa_err]);      
  colormap(flipud(pink(256)))
  axis xy
  colorbar
  title('${\rm Log}_{10}$ of absolute error in $|\nabla u({\bf x})|$', ...
	'Interpreter','LaTeX')
  xlabel('$x_1$','Interpreter','LaTeX')
  ylabel('$x_2$','Interpreter','LaTeX')
  axis(xylim)
  axis equal
%
%
%
  function iout=Ineval(zg,theta,np)
  iout=zeros(np,1);            % point assumed to be inside
  thetag=abs(angle(zg));
  iout(thetag>=theta/2)=1;     % point is outside
  sg=0.5+thetag/theta;
  iout(abs(zg)>=sin(pi*sg))=1; % point is outside

  function [xg,yg,zg,ngrtot]=targinit(xylim,ngr)
  xg=linspace(xylim(1),xylim(2),ngr);
  yg=linspace(xylim(3),xylim(4),ngr);
  ngrtot=ngr^2;
  zg=zeros(ngrtot,1);
  for k=1:ngr
    zg((k-1)*ngr+(1:ngr))=xg(k)+1i*yg;
  end
  
  function [R,Rstor,Kstor]=Rcomp(theta,lambda,LogC,T,W,Pbc,PWbc,...
				 nsub,prema,nse,npan0,iref,starL,circL)
  Rstor=zeros(128,128,nse);
  Kstor=zeros(192,192,nse);
  myind1=[1:48 97:144];
  myind2=[49:96 145:192];
  for level=1:nsub
    [z,zp,nz,kappa,awzp,sidi]=zlocinit(theta,T,W,nsub,level,npan0,iref);  
    MS=MSmatinit(z,nz,kappa,awzp,lambda,96);
    MS(myind1,myind1)=0;
    MS(myind2,myind2)=0;
    ML=MLmatinit(LogC,z,zp,nz,kappa,awzp,sidi,lambda,96);
    MK=MS+ML;
%    
    MAT=eye(192)+MK;
    if level==1                     %  Initializer for R
      R=eye(128);  
    end
    if level >= prema
      Kstor(:,:,level-prema+1)=MK;
      Rstor(:,:,level-prema+1)=R;
    end
    R=SchurBana(Pbc,PWbc,MAT,R,starL,circL);   
  end
%
%
%
  function A=SchurBana(P,PW,K,A,starL,circL)
  starS=[17:48 81:112];
  circS=[1:16 49:80 113:128];
  VA=K(circL,starL)*A;
  PTA=PW'*A;
  PTAU=PTA*K(starL,circL);
  DVAUI=myinv(K(circL,circL)-VA*K(starL,circL));
  DVAUIVAP=DVAUI*(VA*P);
  A(starS,starS)=PTA*P+PTAU*DVAUIVAP;
  A(circS,circS)=DVAUI;
  A(circS,starS)=-DVAUIVAP;
  A(starS,circS)=-PTAU*DVAUI;

  function b=SchurBanaLs(K,A,starL,circL,rhs)
  VA=K(circL,starL)*A;
  AU=A*K(starL,circL);
  DVAU=K(circL,circL)-VA*K(starL,circL);
  b=zeros(length(rhs),1);
  DVAUIrhs=mylins(DVAU,rhs(circL));
  DVAUIVArhs=mylins(DVAU,VA*rhs(starL));
  b(starL)=A*rhs(starL)+AU*DVAUIVArhs-AU*DVAUIrhs;
  b(circL)=DVAUIrhs-DVAUIVArhs;

  function x=mylins(M,b)
% *** equation (9) in Henderson & Searle ***
  np=length(M)/2;
  b1=b(1:np);
  b2=b(np+1:2*np);
  A=M(1:np,1:np);
  U=M(1:np,np+1:2*np);
  V=M(np+1:2*np,1:np);
  D=M(np+1:2*np,np+1:2*np);  
  x=[ (A-U/D*V)\b1-A\(U*((D-V/A*U)\b2));
     -D\(V*((A-U/D*V)\b1))+(D-V/A*U)\b2];

  function Mi=myinv(M)
% *** equation (9) in Henderson & Searle ***
  np=length(M)/2;
  A=M(1:np,1:np);
  U=M(1:np,np+1:2*np);
  V=M(np+1:2*np,1:np);
  D=M(np+1:2*np,np+1:2*np);  
  Mi=[inv(A-U/D*V) -A\U/(D-V/A*U);
      -D\V/(A-U/D*V) inv(D-V/A*U)];
%
%
%
  function [x,it,trueres]=myGMRESR(A,R,b,n,m,tol)
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
%
%
%
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
%
%
%
  function MK=MSmatinit(z,nz,kappa,awzp,lambda,N)
% Construct the smooth part of the full interaction matrix
% for the first Dirichlet
% problem (i.e., u and du/dn are specified on the boundary)
% of the modified biharmonic equation.
  K11=zeros(N);K12=zeros(N);K21=zeros(N);K22=zeros(N);

  nz=nz.';
  for m=1:N
    zm=z(m)-z.';
    
    r1=abs(zm);
    r2=r1.*r1;
    r4=r2.*r2;
    
    rnu=real(conj(nz).*zm);
    rnu2=rnu.*rnu;
    rnu3=rnu.*rnu2;

    rnux=real(conj(nz(m))*zm);
    nunux=real(conj(nz(m))*nz);

    K11(m,:) = rnu3./r4.*(-1+lambda^2*r2/8);
    K12(m,:) = 0.25-0.5*rnu2./r2;

    K21(m,:) = (3*nunux.*rnu2./r2.*(-1+lambda^2*r2/8) ...
		-lambda^2/8*rnux.*rnu ...
		-rnux.*rnu3./r4.*(-4+lambda^2*r2/4-lambda^4*r4/48))./r2;    
    
    K22(m,:) = -(rnu.*nunux + rnux.*rnu2./r2.*(-1+lambda^2*r2/8))./r2;
  end

% diagonal terms
  K11(1:N+1:end) = 0;
  K12(1:N+1:end) = 0.25;
  K21(1:N+1:end) = -3/4*kappa.^2;
  K22(1:N+1:end) = 0.5*kappa;
% multiply by [I 0; 2Kappa I] from right so that 
% the diagonal term is the identity matrix.
  for m=1:N
    K11(:,m)=K11(:,m)+2*kappa(m)*K12(:,m);
    K21(:,m)=K21(:,m)+2*kappa(m)*K22(:,m);
  end
% assemble matrix blocks  
  MK=[K11, K12; K21, K22]; 
% multiply by the weights  
  W=awzp(:,ones(1,N)).';
  W=repmat(W,2,2);
  MK=MK.*W;
  MK=2/pi*MK;
%
%
%
  function MK=MLmatinit(LogC,zs,zp,nz,kappa,awzp,sinterdiff,lambda,N)
% Construct the remaining part of the full interaction matrix for 
% the first Dirichlet problem (i.e., u and du/dn are specified on the boundary)
% of the modified biharmonic equation.
  K11=zeros(N);K12=zeros(N);K21=zeros(N);K22=zeros(N);

  for m=1:N
    zm=zs(m)-zs.';
    
    r1=abs(zm);
    r2=r1.*r1;
    r4=r2.*r2;

    rnu=real(conj(nz.').*zm);
    rnu2=rnu.*rnu;
    rnu3=rnu.*rnu2;

    rnux=real(conj(nz(m))*zm);
    nunux=real(conj(nz(m))*nz.');

    z=lambda*r1;
    % get coefficients c_j (j=0,...,5)
    [c0,c1,c2,c3,c4,c5] = cjseval(z);

    t1=-rnu./r2;
    t2=rnu3./r4;
    K11(m,:) = t1.*c1+t2.*c2;
    
    t0=-1/2+rnu2./r2;
    K12(m,:) = t0.*c0;

    t3=-nunux;
    t4=3*nunux.*rnu2./r2;
    t5=rnux.*rnu./r2;
    t6=-rnux.*rnu3./r4;
    K21(m,:) = (t3.*c1 + t4.*c2 + t5.*c3 + t6.*c4)./r2;    
    
    t7=2*rnu.*nunux;
    t8=rnux/2;
    t9=-rnux.*rnu2./r2;
    K22(m,:) = (t7.*c0 + t8.*c5 + t9.*c2)./r2;
    
    myind = find(LogC(m,:));
    % Now find the logarithmic part of the kernels
    [c0lg,c1lg,c2lg,c3lg,c4lg,c5lg]=cjsLOGeval(z(myind));
    
    K11lg = t1(myind).*c1lg + t2(myind).*c2lg;

    K12lg = t0(myind).*c0lg;
    
    K21lg = ( t3(myind).*c1lg + t4(myind).*c2lg ...
        + t5(myind).*c3lg + t6(myind).*c4lg )./r2(myind);
    
    K22lg = ( t7(myind).*c0lg + t8(myind).*c5lg ...
        + t9(myind).*c2lg )./r2(myind);
    
    % adding the correction for the logarithmic singularity
    K11(m,myind) = K11(m,myind) + K11lg.*LogC(m,myind);
    K12(m,myind) = K12(m,myind) + K12lg.*LogC(m,myind);
    K21(m,myind) = K21(m,myind) + K21lg.*LogC(m,myind);
    K22(m,myind) = K22(m,myind) + K22lg.*LogC(m,myind);
  end
% diagonal terms
  K11(1:N+1:end)=0;
  K12(1:N+1:end)=0;
  K22(1:N+1:end)=0;
  tmp=-lambda^2/8*(log(lambda/2)+0.25+0.5772156649015328606);
  for m=1:N
    dsd=sinterdiff(fix((m-1)/16)+1)/2;
    K21(m,m)=tmp-lambda^2*( LogC(m,m)+log( dsd*abs(zp(m)) ) )/8;
  end
% multiply by [I 0; 2Kappa I] from right so that 
% the diagonal term is the identity matrix.
  for m=1:N
    K11(:,m)=K11(:,m)+2*kappa(m)*K12(:,m);
    K21(:,m)=K21(:,m)+2*kappa(m)*K22(:,m);
  end
% assemble matrix blocks  
  MK=[K11, K12; K21, K22]; 
% multiply by the weights  
  W=awzp(:,ones(1,N)).';
  W=repmat(W,2,2);
  MK=MK.*W;
  MK=2/pi*MK;
%
%
%
  function [c0,c1,c2,c3,c4,c5] = cjseval(z)
  dlim=0.4;
  k0=besselk(0,z);
  k1=besselk(1,z);

  z2=z.^2;
  z4=z2.*z2;
  
  zk1=z.*k1;
  z2k0=z2.*k0;
    
  c0=k0+2./z.*k1-2./z2;
  c1=3*c0+zk1+1/2;
  c2=4*c0+zk1;
  c3=12*c0+5*zk1+z2k0+1;
  c4=24*c0+8*zk1+z2k0;
  c5=2*c0+zk1;

  c0=c0+1/2;
  c2=c2+1-z2/8;
  c3=c3+z2/8;
  c4=c4+4-z2/4+z4/48;
  
  ind = intersect(find(z<dlim),find(z~=0)); % use asymptotic expansion for small z
  if ~isempty(ind)
    gamma=0.57721566490153286061;
    lnz=log(2./z(ind))-gamma;
    z2=z(ind).^2;
    z4=z2.^2;
    z6=z2.*z4;
    z8=z2.*z6;
    z10=z2.*z8;
    z12=z2.*z10;
    z14=z2.*z12;
    
    c0(ind) = z2/8.*((lnz+3/4)+(lnz+17/12).*z2/12+(lnz+43/24).*z4/384 ...
	          +(lnz+247/120).*z6/23040+(lnz+34/15).*z8/2211840 ...
              +(lnz+256/105).*z10/309657600 ...
              +(lnz+1447/560).*z12/59454259200 ...
	          +(lnz+13663/5040).*z14/14982473318400);
    c1(ind) = -z2/8.*((lnz-1/4)+(lnz+13/12).*z2/4 ...
              +(5*lnz+191/24).*z4/384+(7*lnz+1609/120).*z6/23040 ...
              +(lnz+97/45).*z8/245760 ...
              +(11*lnz+2711/105).*z10/309657600 ...
              +(13*lnz+18251/560).*z12/594542549200 ...
              +(lnz+13327/5040).*z14/998831554560);
    c2(ind) = -z4/48.*((lnz+11/12)+(lnz+37/24).*z2/16+(lnz+227/120).*z4/640 ...
              +(lnz+257/120).*z6/46080+(lnz+491/210).*z8/5160960 ...
              +(lnz+4201/1680).*z10/825753600 ...
	          +(lnz+13303/5040).*z12/178362777600 ...
	          +(lnz+2783/1008).*z14/49941577728000);
    c3(ind) = z4/16.*((lnz+7/12)+(5*lnz+161/24).*z2/48 ...
              +(7*lnz+1469/120).*z4/1920+(lnz+731/360).*z6/15360 ...
              +(11*lnz+5191/210).*z8/15482880 ...
              +(13*lnz+52933/1680).*z10/2477260800 ...
              +(lnz+12967/5040).*z12/35672555520 ...
	          +(17*lnz+46303/1008).*z14/149824733184000);
    c4(ind) = z6/384.*((lnz+25/24)+(lnz+197/120).*z2/20+(lnz+79/40).*z4/960 ...
              +(lnz+1859/840).*z6/80640+(lnz+4033/1680).*z8/10321920 ...
	          +(lnz+12883/5040).*z10/1857945600 ...
	          +(lnz+2711/1008).*z12/445906944000 ...
	          +(lnz+31117/11088).*z14/137339338752000);
    c5(ind) = -z2/4.*((lnz+1/4)+1/6*(lnz+7/6).*z2+1/128*(lnz+13/8).*z4 ...
              +(lnz+29/15).*z6/5760+(lnz+13/6).*z8/442368 ...
              +(lnz+989/420).*z10/51609600 ...
              +(lnz+3*67/80).*z12/8493465600 ...
              +(lnz+3337/1260).*z14/1872809164800);
  end
%
%
%
  function [c0lg,c1lg,c2lg,c3lg,c4lg,c5lg]=cjsLOGeval(z)
  dlim=0.4;
  
  k0lg=real(-besselj(0,1i*z)); % log part of besselk(0,z)
  k1lg=real(-1i*besselj(1,1i*z)); % log part of besselk(1,z)
    
  zk1lg=z.*k1lg;
  z2k0lg=z.^2.*k0lg;
    
  c0lg=k0lg+2./z.*k1lg;
  c1lg=3*c0lg+zk1lg;
  c2lg=4*c0lg+zk1lg;
  c3lg=12*c0lg+5*zk1lg+z2k0lg;
  c4lg=24*c0lg+8*zk1lg+z2k0lg;
  c5lg=2*c0lg+zk1lg;
    
  ind = intersect(find(z<dlim),find(z~=0)); % use Taylor expansion for small z
  if ~isempty(ind)
    z1 = z(ind);
    z2 = z1.^2;
    z4 = z2.^2;
    z6 = z2.*z4;
    z8 = z2.*z6;
    z10=z2.*z8;
    z12=z2.*z10;
    z14=z2.*z12;
    
    c0lg(ind) = -z2/8.*(1+z2/12+z4/384+z6/23040+z8/2211840 ...
                +z10/309657600 + z12/59454259200 ...
                +z14/14982473318400); 
    c1lg(ind) = z2/8.*(1+z2/4+5*z4/384+7*z6/23040+z8/245760 ...
                +11*z10/309657600 + 13*z12/59454259200 ...
                +z14/998831554560); 
    c2lg(ind) = z4/48.*(1+z2/16+z4/640+z6/46080+z8/5160960 ...
                +z10/825753600 + z12/178362777600 ...
                +z14/49941577728000); 
    c3lg(ind) = -z4/16.*(1+5*z2/48+7*z4/1920+z6/15360+11*z8/15482880 ...
                +13*z10/2477260800 + z12/35672555520 ...
                +17*z14/149824733184000);
    c4lg(ind) = -z6/384.*(1+z2/20+z4/960+z6/80640+z8/10321920 ...
                + z10/1857945600 + z12/445906944000 ...
                + z14/137339338752000); 
    c5lg(ind) = z2/4.*(1+z2/6+z4/128+z6/5760+z8/442368 ...
                +z10/51609600 + z12/8493465600 ...
                +z14/1872809164800);
  end
%
%
%   
  function [z,zp,nz,wzp,kappa,awzp,zinter,np,s]=zinit(theta,sinter, ...
						    sinterdiff,T,W,npan)
  np=16*npan;
  s=zeros(np,1);
  w=zeros(np,1);
  for k=1:npan
    myind=(k-1)*16+1:k*16;
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
  nz(np/2+1:np)=conj(flipud(nz(1:np/2)));  
  wzp=w.*zp;
  awzp=w.*abs(zp);
  kappa = imag(zpp./zp)./abs(zp);

  function [z,zp,nz,kappa,awzp,sidi]=zlocinit(theta,T,W,nsub,level,npan0,iref)
  denom=2^(nsub-level)*npan0*2^iref;
  s=[T/4+0.25;T/4+0.75;T/2+1.5]/denom;
  w=[W/4;W/4;W/2]/denom;
  w=[flipud(w);w];
  sidi=[1;0.5;0.5;0.5;0.5;1]/denom;
  z=zfunc(s,theta);
  z=[conj(flipud(z));z];
  zp=zpfunc(s,theta);
  nz=-1i*zp./abs(zp);
  zp=[-conj(flipud(zp));zp];
  nz=[conj(flipud(nz));nz];
  zpp=zppfunc(s,theta);
  zpp=[conj(flipud(zpp));zpp];
  
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
%
%
%
  function M1=LogCinit(sinterdiff,T,W,npan,N,iper)
% *** excluding diagonal entries ***
% *** Corrections to Logarithmic potential log(|tau-z|) ***   
% *** block-tri-diagonal output ***
% *** iper=0,1 (0 is open arc, 1 is closed contour) ***
  M1=zeros(N);
  A=ones(16);
  for k=2:16
    A(:,k)=T.*A(:,k-1);
  end
  if iper==1
    kstart=1;
    kend=npan;  
  else
    kstart=2;
    kend=npan-1;  
  end
% *** central blocks ***
  TMP=-log(abs(T(:,ones(1,16))-T(:,ones(1,16))'));
  TMP(1:17:256)=0;
  TMP=TMP+WfrakLinit(A,0,1,T)./W(:,ones(1,16))';
  for k=1:npan
    myind=(k-1)*16+1:k*16;
    M1(myind,myind)=TMP;
  end
% *** superdiagonal blocks (targets to the left) ***
  for k=kstart:npan
    myinds=(k-1)*16+1:k*16;
    km1=mod(k-2,npan)+1;
    mi=(km1-1)*16+1:km1*16;       
    alpha=sinterdiff(km1)/sinterdiff(k);          
    [TMP,accept,na]=WfrakLinit(A,-1-alpha,alpha,T);
    mi=mi(accept);
    for nj=1:16
      M1(mi,myinds(nj))=-log(abs(T(nj)+1+(1-T(accept))*alpha));       
    end
    M1(mi,myinds)=M1(mi,myinds)+TMP./W(:,ones(1,na))';
  end
% *** subdiagonal blocks (targets to the right) ***
  for k=1:kend
    myinds=(k-1)*16+1:k*16;
    kp1=mod(k,npan)+1;
    mi=(kp1-1)*16+1:kp1*16;       
    alpha=sinterdiff(kp1)/sinterdiff(k);               
    [TMP,accept,na]=WfrakLinit(A,1+alpha,alpha,T);
    mi=mi(accept);
    for nj=1:16
      M1(mi,myinds(nj))=-log(abs(T(nj)-1-(T(accept)+1)*alpha));   
    end
    M1(mi,myinds)=M1(mi,myinds)+TMP./W(:,ones(1,na))';
  end
  if iper==1
    M1=sparse(M1);
  end
%
%
%
  function [WfrakL,accept,na]=WfrakLinit(A,trans,mscale,T)
% *** T is target vector, sources on canonical interval ***
  accept=1:16;
  T=trans+mscale*T;
  accept=accept(abs(T)<2);
  na=length(accept);
  p=zeros(1,17);
  Q=zeros(16,na);
  c=(1-(-1).^(1:16))./(1:16);
  for j=1:na
    jj=accept(j);
    p(1)=log(abs((1-T(jj))/(1+T(jj))));
    p111=log(abs(1-T(jj)^2));
    for k=1:16
      p(k+1)=T(jj)*p(k)+c(k);
    end
    Q(1:2:15,j)=(p111-p(2:2:16))./(1:2:15);
    Q(2:2:16,j)=(p(1)-p(3:2:17))./(2:2:16);
  end
  WfrakL=Q.'/A;
%
%
%
  function [IP,IPW]=IPinit(T,W)
  A=ones(16);
  AA=ones(32,16);
  T2=[T-1;T+1]/2;
  W2=[W;W]/2;
  for k=2:16
    A(:,k)=A(:,k-1).*T;
    AA(:,k)=AA(:,k-1).*T2;   
  end
  IP=AA/A;
  IPW=IP.*(W2*(1./W)');
%
%
%
  function Pbc=Pbcinit(IP)
  Pbc=blkdiag(IP,rot90(IP,2));
%
%
%
  function Pbc=Poper(IP)
  Pbc=blkdiag(eye(16),IP,rot90(IP,2),eye(16));
%
%
%
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

  function [f,fz]=testfun(ztarg,lambda)
% *** the fundamental solution of the modified biharmonic equation in 2D ***
  persistent zs ns q
  q=[0.171;0.720; 0.918; 0.825; 0.078];
  x=[1.280;0.413;-0.235;-0.145; 1.050];
  y=[0.177;0.795; 0.315;-0.474;-0.580];
  zs=x+1i*y;
  ns=5;
%
  f = zeros(size(ztarg));
  if nargout>1 
    fz = zeros(size(ztarg)); 
  end
%    
  for i=1:ns
    zi=ztarg-zs(i);
    r1=abs(zi);
    z=lambda*r1;
    f = f+q(i)*(log(r1)+besselk(0,z)); 
    if nargout>1
      fz =fz+q(i)*(zi./r1.^2-lambda*zi.*besselk(1,z)./r1);
    end
  end

  function [wcmpL,wcmpC,wcmpH,wcmpS]=wLCHSinit(a,b,ztg,zsc,nz,wzp,awzp,io)
% *** zt is target vector in transformed plane ***
  Ng=length(zsc);
  cc=(b-a)/2;
  zt=(ztg-(b+a)/2)/cc;
  zsctr=(zsc-(b+a)/2)/cc;
%  c1=(1-(-1).^(1:Ng))./(1:Ng);
  c2=-(1/(1-zt)+(-1).^(0:Ng-1)/(1+zt));
  c3=-0.5*(1/(1-zt)^2+(-1).^(1:Ng)/(1+zt)^2)...
    -0.25*(0:Ng-1).*(1:Ng).*(1/(1-zt)+(-1).^(1:Ng)/(1+zt));
  upp=log(1-zt);
  loo=log(-1-zt);
  if io==1 % point is outside
    if imag(zt)>0 && abs(real(zt))<1
      loo=loo+2i*pi;
    end
  else     % point is inside
    if imag(zt)<0 && abs(real(zt))<1
      loo=loo-2i*pi;    
    end
  end
  f=zeros(Ng+1,1); % f(k+1)=\int_{-1}^1 p_k(z)/(z-zt) dz 
  fp=zeros(Ng,1);  % fp(k+1)=\int_{-1}^1 p'_k(z)/(z-zt) dz 
  fpp=zeros(Ng,1); % fpp(k+1)=\int_{-1}^1 p"_k(z)/(z-zt) dz 
  q=zeros(Ng,1);   % q(k+1)=\int_{-1}^1 p_k(z)*log(z-zt) dz 
  f(1)=upp-loo;
  f(2)=2+zt*f(1);
  fp(2)=f(1);
  for k=1:Ng-1
    f(k+2)  =( (2*k+1)*zt*f(k+1)-k*f(k) )/(k+1);
  end
  q(1)=2*upp-(zt+1)*f(1)-2;
  for k=1:Ng-1
    q(k+1)=(f(k)-f(k+2))/(2*k+1);
  end

  for k=1:Ng-2
    fp(k+2) =(2*k+1)*f(k+1)+fp(k);
    fpp(k+2)=(2*k+1)*fp(k+1)+fpp(k);
  end
% Now add corrections so that fp(k+1)=\int_{-1}^1 p_k(z)/(z-zt)^2 dz
% with p_k the Legendre polynomial of degree k.
  fp=fp+c2.'; 
% Now add corrections so that fpp(k+1)=\int_{-1}^1 p_k(z)/(z-zt)^3 dz
  fpp=0.5*fpp+c3.';
  U=legepols(zsctr,Ng-1);
  U=inv(U.');
  wcmpL=imag(U*q*cc.*conj(nz))./awzp-log(abs((zsc-ztg)/cc));
  wcmpC=U*f(1:Ng)/1i-wzp./(zsc-ztg)/1i;
  wcmpH=U*fp/1i/cc-wzp./(zsc-ztg).^2/1i;
  wcmpS=U*fpp/1i/cc.^2-wzp./(zsc-ztg).^3/1i;
%   U=U.';
%   wcmpL=imag(U\q*cc.*conj(nz))./awzp-log(abs((zsc-ztg)/cc));
%   wcmpC=U\f(1:Ng)/1i-wzp./(zsc-ztg)/1i;
%   wcmpH=U\fp/1i/cc-wzp./(zsc-ztg).^2/1i;
%   wcmpS=U\fpp/1i/cc.^2-wzp./(zsc-ztg).^3/1i;

  function p = legepols(x,n)
  % return values of legendre polynomial p_0(x), p_1(x), ..., p_n(x)
  m=length(x);
  p=ones(m,n+1);
  if n == 0 
    return
  end
  p(:,2)=x;
  for k=1:n-1
    p(:,k+2)=( (2*k+1)*x.*p(:,k+1) - k*p(:,k) )/(k+1);
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
%while niter < nmax && (err>eps || abs(w3)>1e-30),
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

  function f = mullerfunc(x,theta,zt)
% return the value of the Legendre expansion of order n at x,
% the expansion coefficients are given by coef.
  z=zfunc(x,theta);
  zp=zpfunc(x,theta);
  zn=-1i*zp./abs(zp);
  f=imag((zt-z)*conj(zn));

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
  