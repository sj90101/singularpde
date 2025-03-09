  function modgp1j
% *****************************************************************
% * Produces Figure 11 in                                         *
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
% asymp: set asymp=1 if asymptotics is of prime concern
%        set asymp=0 if field evaluation is of prime concern
%
% iref:  Coarse mesh refinement towards corner vertex. Recommended
%        value is iref=0. A larger iref, in combination with dlim1=1,
%        could lead to speedups in the post-processor.
%
% *** User specified quantities ********************
  lambda=20  % parameter in the modified biharmonic equation
  nv=3;      % number of vertices
  theta=[2.58e-01 2.20 4.27];
  zv=exp(1i*theta'); % coordinates of vertices in complex number form
  lmin=min(abs(wshift('1D',zv,1)-zv));
  npan0min=24;
  hc=lmin/npan0min; % length of two nearest panels to every corner
  nsub=40        % number of subdivisions
  asymp=0;       % see explanation above
  dlim1=0.8      % closeness parameter for 16-point quadrature
  ngr=300;       % field evaluation at ngr^2 points
  xylim=[-1 1.1 -1 1];  % computational domain
% **************************************************
  Ng=16;
  T=Tinit16;
  W=Winit16;
  [IP,IPW]=IPinit(T,W);
  Pbc=Pbcinit(IP);
  PWbc=Pbcinit(IPW);
  Pbc2=blkdiag(Pbc,Pbc);
  PWbc2=blkdiag(PWbc,PWbc);
%
  sidiloc=[1;0.5;0.5;0.5;0.5;1];
  LogCloc=LogCinit(sidiloc,T,W,6,6*Ng,0);
%
% *** Panels, points, and weights ***
  [sinter,sinterdiff,npan,sc]=panelinit(hc,zv,nv);
  [z,zp,nz,wzp,kappa,awzp,zinter,np]=zinit(sinter,sinterdiff,T,W,npan,zv,nv);
% total number of discretization points
  nps=sum(np);
  ntot=2*nps; % total number of unknowns/equations
  [xg,yg,zg,ngrtot]=targinit(xylim,ngr);
%
% *** The whole interaction matrix M_coa^{\circ} is set up ***
  sinterdiffs=[];
  for i=1:nv
    sinterdiffs=[sinterdiffs; sinterdiff{i}];
  end
  npans=sum(npan);
  LogC=LogCinit(sinterdiffs,T,W,npans,nps,1);
% *** smooth part of MAT
  MS=MSmatinit(z,nz,kappa,awzp,lambda,nps); 
% *** log part of MAT
  ML=MLmatinit(LogC,z,zp,nz,kappa,awzp,sinterdiffs,lambda,nps,Ng); 
%
% set all nonsmooth blocks of MS to zero
%
% vertex #1  
  bl=-2*Ng+1:0;
  br=1:2*Ng;
  myind1=[nps+bl 2*nps+bl ];
  myind2=[br nps+br];
  MS(myind1,myind2)=zeros(4*Ng);
  MS(myind2,myind1)=zeros(4*Ng);  
% rest vertices
  for iv=2:nv
    npi=sum(np(1:iv-1));
    myind1=[npi+bl nps+npi+bl];
    myind2=[npi+br nps+npi+br];
    MS(myind1,myind2)=zeros(4*Ng);
    MS(myind2,myind1)=zeros(4*Ng);  
  end
% set ML^\star to zero
  starind=zeros(nv,8*Ng);
  bt=-2*Ng+1:2*Ng;
  for iv=1:nv
    if iv>1
      npi=sum(np(1:iv-1));
      starind(iv,:)=[npi+bt nps+npi+bt];
    else
      starind1=[nps-2*Ng+1:nps 1:2*Ng];
      starind1=[starind1 nps+starind1];
      starind(iv,:)=starind1;
    end
  end    
  for iv=1:nv
    ML(starind(iv,:),starind(iv,:))=zeros(8*Ng);  
  end
  Mcirc=MS+ML;  
%
% *** Construct the right hand side ***
  [g,gz]=testfun(z,lambda);
  gn=real(conj(nz).*gz);
  rhs=[g;gn];  
%  
% *** Initializer ***
  starL=[Ng+1:5*Ng 7*Ng+1:11*Ng];
  circL=[1:Ng 5*Ng+1:7*Ng 11*Ng+1:12*Ng];
  Rc=zeros(8*Ng,8*Ng,nv);
  t1=cputime;
  for iv=1:nv
      Rc(:,:,iv)=initializeRc(T,W,zv,nv,iv,Pbc2,PWbc2,starL,circL);
  end
  t2=cputime-t1;
  disp(['The CPU time on the initializer is ',num2str(t2), ' seconds'])
%
  if asymp==1
    prema=2
  else
    prema=min(nsub-ceil(log2(sum(awzp(1:Ng))/min(min(abs(zg-zv.'))))),nsub)
  end
  disp(['nsub = ',num2str(nsub)])
  %prema=min(nsub-ceil(log2(sum(awzp(1:16))/min(abs(zg)))),nsub);
  nse=nsub-prema+1;
  disp(['nse = ',num2str(nse)])
  t1=cputime;
  RG=speye(ntot);
  Rstor=zeros(8*Ng,8*Ng,nse,nv);
  Kstor=zeros(12*Ng,12*Ng,nse,nv);
  for iv=1:nv
% *** recursion for the R matrix for each vertex***
    [R,Rstor(:,:,:,iv),Kstor(:,:,:,iv)]=Rcomp(Rc(:,:,iv),zv,nv,iv, ...
		lambda,LogCloc,T,W,Pbc2,PWbc2,nsub,prema,sc,starL,circL);
    RG(starind(iv,:),starind(iv,:))=R;  
  end
  t2=cputime-t1;
  disp(['The CPU time on the recursion is ',num2str(t2), ' seconds'])
  t1=cputime;
  [rhotilde,iter]=myGMRESR(Mcirc,RG,rhs,ntot,100,5*eps);
  t2=cputime-t1;
  disp(['The CPU time on GMRES is ',num2str(t2), ' seconds'])
  disp(['GMRES iter RCIP = ',num2str(iter)]) 
%
% *** Construction of fine grid ***
  [sinterfin,sinterdifffin,npanfin]=cornerrefine(sinter,sinterdiff,npan,nse,nv);
  [zfin,~,nzfin,wzpfin,kappafin,awzpfin,zinterfin,npfin,sfin]= ...
      zinit(sinterfin,sinterdifffin,T,W,npanfin,zv,nv);
%
% *** Reconstruction of layer density on the fine grid ***
  npsfin=sum(npfin);
  rhofin=zeros(2*npsfin,1);
  for iv=1:nv
    rhot=rhotilde(starind(iv,:));
    pts2=Ng*(nse+2);
    rhofinloc=zeros(4*pts2,1);
    PbcL=Poper(IP);
    Pbc2L=blkdiag(PbcL,PbcL);
    for k=nse:-1:1
      Kcirc=Kstor(:,:,k,iv);
      Kcirc(starL,starL)=zeros(8*Ng);
      MAT=eye(12*Ng)+Kcirc;
      tmp=Pbc2L*rhot;
      rhot=tmp-Kcirc*SchurBanaLs(MAT,Rstor(:,:,k,iv),starL,circL,tmp,nse-k+1);
      myindL=[pts2-(k+2)*Ng+(1:Ng) pts2+(k+1)*Ng+(1:Ng)];
      rhofinloc([myindL 2*pts2+myindL])=rhot(circL);
      rhot=rhot(starL);
    end
    myindL=nse*Ng+1:(4+nse)*Ng;
    rhofinloc([myindL 2*pts2+myindL])=Rstor(:,:,1,iv)*rhot;
% now put everything in the right place of rhofin
    ind=1:np(iv)-4*Ng; % points on the middle panels
    if iv == 1
      ind1=1:pts2;
      ind2=-pts2+1:0;
% density 1    
      rhofin(ind1)=rhofinloc(pts2+ind1);
% middle panels are unchanged, so simply copy from rhotilde
      rhofin(pts2+ind)=rhotilde(2*Ng+ind); 
      rhofin(npsfin+ind2)=rhofinloc(ind1);
% density 2    
      rhofin(npsfin+ind1)=rhofinloc(3*pts2+ind1);
% middle panels are unchanged, so simply copy from rhotilde
      rhofin(npsfin+pts2+ind)=rhotilde(nps+2*Ng+ind); 
      rhofin(2*npsfin+ind2)=rhofinloc(2*pts2+ind1);
    else
      indcloc=1:2*pts2;
      indcfin=-pts2+1:pts2;
      ic=sum(np(1:iv-1));
      icfin=sum(npfin(1:iv-1));
% density 1
      rhofin(icfin+indcfin)=rhofinloc(indcloc);
      rhofin(icfin+pts2+ind)=rhotilde(ic+2*Ng+ind);
% density 2
      rhofin(npsfin+icfin+indcfin)=rhofinloc(2*pts2+indcloc);
      rhofin(npsfin+icfin+pts2+ind)=rhotilde(nps+ic+2*Ng+ind);
    end      
  end
  rho1fin=rhofin(1:npsfin);
  rho2fin=rhofin(npsfin+1:2*npsfin);
%
% *** Zoom of layer densities close to corner ***
  if asymp==1
    figure(1)
    iv=1;
    if iv==1
      npi=0;
    else
      npi=sum(npfin(1:iv-1));
    end
% *** note that points 1:32 must be omitted. They are not rho, but rhohat
    myind=2*Ng+1:npfin(iv)/2;
    loglog(sfin{iv}(myind),abs(rho1fin(npi+myind)),'b.', ...
	   sfin{iv}(myind),abs(rho2fin(npi+myind)),'r.')
    title('Reconstruction of $\rho$ to the right of the corner', ...
	  'Interpreter','LaTeX')
    grid on
    xlabel('distance to corner','Interpreter','LaTeX')
    hl=legend('$|\rho_1|$','$|\rho_2|$','Location','NorthEast');
    set(hl,'Interpreter','LaTeX','FontSize',12)
    drawnow
  end
%
% *** Evaluation of the field and its gradient ***
  rhohat=RG*rhotilde;
  rho1 = rhohat(1:nps);
  rho2 = rhohat(nps+1:end);
  [uref,uzref]=testfun(zg,lambda); % reference solution
  iout=Ineval(zg,zv,ngrtot);
  t1=cputime;
  [u,uz]=fieldcomp(lambda,rho1,rho2,z,nz,wzp,awzp,kappa,zinter,rho1fin, ...
	rho2fin,zfin,nzfin,wzpfin,awzpfin,kappafin,zinterfin,sinterfin,...
	zg,dlim1,ngrtot,iout,npan,npanfin,Ng,zv,nv);
  t2=cputime-t1;
  disp(['The CPU time on evaluation is ',num2str(t2), ' seconds'])
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
    figure(6)
    plot(real(z),imag(z),'b',real(zg(ind)),imag(zg(ind)),'r*',real(zinter),imag(zinter),'g.')
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

  function [u,uz]=fieldcomp(lambda,rho1,rho2,z,nz,wzp,awzp,kappa,zinter, ...
	rho1fin,rho2fin,zfin,nzfin,wzpfin,awzpfin,kappafin,zinterfin, ...
	sinterfin,zg,dlim1,ngrtot,iout,npan,npanfin,Ng,zv,nv)
  dlim4=0.04;
% *** Initial screening ***
  disp('Screening starts')
% *** 1 => regular evaluation
% *** 2 => some close evaluation
% *** 3 => some close corner-evaluation
  npans=sum(npan); % total number of panels
  cpi=[1 npans]; % corner panel indices
  for iv=2:nv
    npi=sum(npan(1:iv-1));
    cpi=[cpi npi npi+1];
  end
  panlen=zeros(npans,1);
  for k=1:npans
    panlen(k)=sum(awzp((k-1)*Ng+1:k*Ng));
  end
  npansfin=sum(npanfin);
  panlenfin=zeros(npansfin,1);
  for k=1:npansfin
    panlenfin(k)=sum(awzpfin((k-1)*Ng+1:k*Ng));
  end
  kvec=find(iout==0)';          % indicies of points inside
  mylist=zeros(ngrtot,1); 
  mylist(kvec)=1;
  for k=kvec
    for kk=cpi
      myind=(kk-1)*Ng+1:kk*Ng;
      d=min(abs(z(myind)-zg(k)))/panlen(kk);
      if d<2.2
	mylist(k)=mylist(k)+2;  % close corner-evaluation needed
	break
      end
    end     
    for kk=1:npans
      myind=(kk-1)*Ng+1:kk*Ng;
      d=min(abs(z(myind)-zg(k)))/panlen(kk);
      if d<dlim1 && mylist(k)==1
	mylist(k)=mylist(k)+1;  % special evaluation needed
	break
      end
    end
    for kk=1:npansfin
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
  [zclose,dclose]=findbdpoint(zg(mylist4),sinterfin,npanfin,panlenfin,zv,nv);
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
		   dlim1,panlen,npan,0,Ng);
  end
  disp(['List3 starts, target points = ',num2str(length(mylist3))])
  for k=mylist3
    [u(k),uz(k)]=speceval(lambda,zg(k),zfin,nzfin,wzpfin,kappafin, ...
	awzpfin,zinterfin,rho1fin,rho2fin,dlim1,panlenfin,npanfin,0,Ng);
  end
  disp(['List4 starts, target points = ',num2str(length(mylist4))])
  k4=1;
  for k=mylist4
    [u(k),uz(k)]=closeeval(lambda,zg(k),zclose(k4),dclose(k4),zfin, ...
	nzfin,wzpfin,kappafin,awzpfin,zinterfin,rho1fin,rho2fin,dlim1, ...
	panlenfin,npanfin,0,Ng);
    k4=k4+1;
  end 
  u(iout==1)=mean(u(kvec)); % outside values set to something

  function [zc,dc]=findbdpoint(zt,sinter,npan,panlen,zv,nv)
% find the closed point zc, the shortest distance relative to the 
% panel length, the corresponding panel index and side index to the 
% given target points zt.
  nt=length(zt);
  zc=zeros(nt,1);
  dc=zeros(nt,1);
  kc=zeros(nt,1);
  iseg=zeros(nt,1);
  dz=diff([zv;zv(1)]);
  for k=1:nt
    sp=real((zt(k)-zv).*conj(dz))./abs(dz).^2;
    zcs=zv+dz.*sp; % the point closest to zt(k) on each side
    ds=abs(zt(k)-zcs); % the shortest distance from zt(k) to each side
    dsrel=ds; 
    pani=zeros(nv,1);
    for i=1:nv
      for j=1:npan(i)
	if sp(i)>sinter{i}(j) && sp(i)<=sinter{i}(j+1)
	  pani(i)=j; % the panel index on side i where the closest point lies. 
	  break
	end
      end
      if i==1
	si=0;
      else
	si=sum(npan(1:i-1));
      end
      dsrel(i)=dsrel(i)/panlen(si+pani(i)); 
    end
    [~,id]=min(dsrel);
    zc(k)=zcs(id);
    dc(k)=dsrel(id);
    kc(k)=pani(id);
    iseg(k)=id;
  end    

  function [u,uz]=regueval(lambda,zt,z,nz,kappa,awzp,rho1,rho2)
  [K1,K2,K1z,K2z]=repkernels(zt,z,nz,kappa,lambda);
  u=awzp'*(rho1.*K1+rho2.*K2);
  u=u*2/pi;

  uz=awzp'*(rho1.*K1z+rho2.*K2z);
  uz=uz*2/pi;

  function [u,uz]=speceval(lambda,zt,z,nz,wzp,kappa,awzp,zinter,rho1,rho2,...
                       dlim1,panlen,npan,io,Ng)
  [u,uz]=regueval(lambda,zt,z,nz,kappa,awzp,rho1,rho2);
  npans=sum(npan);
  for k=1:npans
    myind=(k-1)*Ng+(1:Ng);
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


  function [u,uz]=closeeval(lambda,zt,zclose,dclose,z,nz,wzp,kappa,awzp, ...
			    zinter,rho1,rho2,dlim1, panlen,npan,io,Ng)
  n=8;
  x=(cos(pi*(n:-1:0)'/n)+1)/2*0.8; % Chebyshev nodes in ascending order  
  zts=zclose+x*(zt-zclose)/dclose;
  f=zeros(n+1,1);
  [~,f(1)]=testfun(zclose,lambda); % value at the closest boundary point
  for k=2:n+1
    [~,f(k)]=speceval(lambda,zts(k),z,nz,wzp,kappa,awzp,zinter, ...
                      rho1,rho2,dlim1,panlen,npan,io,Ng);
  end
  uz=chebbary(dclose,x,f);
  [u,~]=speceval(lambda,zt,z,nz,wzp,kappa,awzp,zinter,rho1,rho2,dlim1, ...
 			   panlen,npan,io,Ng);

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
%
%
%
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
  F1=zeros(ngr);
  for k=1:ngr
    F1(1:ngr,k)=u((k-1)*ngr+(1:ngr));
  end
  fh =  findobj('type','figure');
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
  fh =  findobj('type','figure');
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
  F1=zeros(ngr);
  for k=1:ngr
    F1(1:ngr,k)=abs(uref((k-1)*ngr+(1:ngr))-u((k-1)*ngr+(1:ngr)));
  end
  F1=log10(F1);
  maxa_err=max(max(F1));
  disp(['max field error = ',num2str(max(max(F1)))])
  fh =  findobj('type','figure');
  nf = length(fh);
  set(groot,'DefaultAxesFontSize',15)
  figure(nf+1)
  imagesc(xg,yg,F1,[log10(eps) maxa_err]);      
  colormap(flipud(pink(256)))
  axis xy
  colorbar
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
  xlabel('$x_1$','Interpreter','LaTeX')
  ylabel('$x_2$','Interpreter','LaTeX')
  axis(xylim)
  axis equal
%
%
%
  function iout=Ineval(zg,zv,np)
  iout=ones(np,1);            % point assumed to be outside

  xv=real([zv; zv(1)]);
  yv=imag([zv; zv(1)]);
  xg=real(zg);
  yg=imag(zg);
  
  [in,on]=inpolygon(real(zg),imag(zg),xv,yv);
  iout(in)=0;

  function [xg,yg,zg,ngrtot]=targinit(xylim,ngr)
  xg=linspace(xylim(1),xylim(2),ngr);
  yg=linspace(xylim(3),xylim(4),ngr);
  ngrtot=ngr^2;
  zg=zeros(ngrtot,1);
  for k=1:ngr
    zg((k-1)*ngr+(1:ngr))=xg(k)+1i*yg;
  end

  function [R,Rstor,Kstor]=Rcomp(Rc,zv,nv,iv,lambda,LogC,T,W,Pbc,PWbc, ...
	nsub,prema,sc,starL,circL)
  Ng=length(T);
  myind1=[1:3*Ng 6*Ng+1:9*Ng];
  myind2=[3*Ng+1:6*Ng 9*Ng+1:12*Ng];
  denom=2^(nsub-1)/sc(iv);
  D21=[ones(4*Ng) ones(4*Ng)/denom;ones(4*Ng)*denom ones(4*Ng)];
  R=D21.*Rc;
% *** Recursion ***
  for level=1:nsub
    [z,zp,nz,~,kappa,awzp,sidi]=zlocinit(zv,nv,iv,T,W,nsub,level,sc);  
    MS=MSmatinit(z,nz,kappa,awzp,lambda,6*Ng);
    MS(myind1,myind1)=0;
    MS(myind2,myind2)=0;
    ML=MLmatinit(LogC,z,zp,nz,kappa,awzp,sidi,lambda,6*Ng,Ng);
    MK=MS+ML;
    MAT=eye(12*Ng)+MK;
    if level >= prema
      Kstor(:,:,level-prema+1)=MK;
      Rstor(:,:,level-prema+1)=R;
    end
    
    R=SchurBana(Pbc,PWbc,MAT,R,starL,circL,nsub-level+1);   
  end

  function R=initializeRc(T,W,zv,nv,iv,Pbc,PWbc,starL,circL)
  Ng=length(T);
  myind1=[1:3*Ng 6*Ng+1:9*Ng];
  myind2=[3*Ng+1:6*Ng 9*Ng+1:12*Ng];
  R=eye(8*Ng);
  myerr=1;
  iter=0;
  D21=[ones(4*Ng) 0.5*ones(4*Ng);2*ones(4*Ng) ones(4*Ng)];
  [z,~,nz,wzp,kappa,awzp,~]=zlocinit(zv,nv,iv,T,W,1,1,ones(nv,1));  
  MK=MKmatinit(z,nz,wzp,kappa,awzp,6*Ng);
  MK(myind1,myind1)=0;
  MK(myind2,myind2)=0;
  MAT=eye(12*Ng)+MK;
  while myerr>2*eps && iter<2000
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
  Mi=d2(:,ones(1,2*np)).*((d1*d2'.*M)\diag(d1));

  function x=mylinsnew(M,b,idepth)
% *** simple diagonal scaling ***
  np=length(M)/2;
  d1=[ones(np,1);ones(np,1)/2^idepth];
  d2=[ones(np,1);ones(np,1)*2^idepth];
  x=d2.*((d1*d2'.*M)\(d1.*b));  
  
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
  function MK=MLmatinit(LogC,zs,zp,nz,kappa,awzp,sinterdiff,lambda,N,Ng)
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
    dsd=sinterdiff(fix((m-1)/Ng)+1)/2;
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

  function [z,zp,nz,wzp,kappa,awzp,zinter,np,s]=zinit(sinter,sinterdiff, ...
						      T,W,npan,zv,nv)
  Ng=length(T);
  np=zeros(nv,1);
  s=cell(nv,1);
  z=[];
  zp=[];
  zpp=[];
  w=[];
  zinter=[];
  for i=1:nv
    np(i)=Ng*npan(i);
    s{i}=zeros(np(i),1);
    wi=zeros(np(i),1);
    for k=1:npan(i)
      myind=(k-1)*Ng+1:k*Ng;
      sdif=sinterdiff{i}(k)/2;
      s{i}(myind)=(sinter{i}(k)+sinter{i}(k+1))/2+sdif*T;
      wi(myind)=W*sdif;
    end
    [zi,zpi,zppi]=zfunc(s{i},zv,nv,i);
    z=[z; zi];
    zp=[zp; zpi];
    zpp=[zpp; zppi];
    w=[w; wi];
    [zinteri,~,~]=zfunc(sinter{i},zv,nv,i);
    if i<nv
      zinteri=zinteri(1:end-1);
    end
    zinter=[zinter; zinteri];
  end
% *************************************************
  nz=-1i*zp./abs(zp);
  wzp=w.*zp;
  awzp=w.*abs(zp);
  kappa = imag(zpp./zp)./abs(zp);
 
  function [z,zp,nz,wzp,kappa,awzp,sidi]=zlocinit(zv,nv,iv,T,W,nsub,level,sc)
% points around the iv-th vertex in counterclockwise order
    denom=2^(nsub-level);
    s=[T/4+0.25;T/4+0.75;T/2+1.5]/denom*sc(iv);
    w=[W/4;W/4;W/2]/denom*sc(iv);
    w=[flipud(w);w];
    sidi=[1;0.5;0.5;0.5;0.5;1]/denom*sc(iv);

    v1=zv(iv);
    if iv==nv
      v2=zv(1);
    else
      v2=zv(iv+1);
    end    
    
    if iv==1
      v3=zv(nv);
    else
      v3=zv(iv-1);
    end

    theta=angle((v3-v1)/(v2-v1))/2;
    
    z=abs((v2-v1))*s*exp(-1i*theta);
    zp=abs((v2-v1))*exp(-1i*theta)*ones(size(z));
    
    z=[conj(flipud(z));z];
    zp=[-conj(flipud(zp));zp];
    
    nz=-1i*zp./abs(zp);
    wzp=w.*zp;
    awzp=w.*abs(zp);
    kappa = zeros(size(zp));
    
  function [z,zp,zpp]=zfunc(s,zv,nv,iseg)
  % zv: coordinates of vertices in complex number form
  % nv: number of vertices   
  % iseg: the ith segment of the polygon
    z1=zv(iseg);
    if iseg<nv
      z2=zv(iseg+1);
    else
      z2=zv(1);
    end
  % change here to have a general curve for each piece.
  % the parametrization is always assumed to be from 0 to 1 for every piece.   
    z=z1+(z2-z1)*s;
    zp=(z2-z1)*ones(size(s));
    zpp=zeros(size(s));
%
%
%
  function [sinter,sinterdiff,npanfin,sc]=panelinit(hc,zv,nv)
    sinter=cell(nv,1);
    sinterdiff=cell(nv,1);
    npanfin=zeros(nv,1);
    sc=zeros(nv,1);
    
    for i=1:nv
      [sinter{i},sinterdiff{i},npanfin(i),sc(i)]=panelinit0(hc,zv,nv,i);
    end
    
  function [sinter,sinterdiff,npanfin,sc]=panelinit0(hc,zv,nv,iseg)
    z1=zv(iseg);
    if iseg==nv
      z2=zv(1);
    else
      z2=zv(iseg+1);
    end

    l=abs(z2-z1);
    nmid=round((l-4*hc)/hc);
    sc=hc/l;
    
    npanfin=nmid+4;
    sinter=zeros(npanfin+1,1);

    sinter(1:3)=(0:2)*sc; % the first two panels
    sinter(npanfin-1:npanfin+1)=1-(2:-1:0)*sc; % the last two panels
    smid=(1-4*sc)/nmid;
    sinter(4:npanfin-1)=2*sc+(1:nmid)*smid; % mid panels

    sinterdiff=diff(sinter);
%
%
%
  function [sinterfin,sinterdifffin,npanfin]=cornerrefine(sinter,sinterdiff, ...
							  npan,nsub,nv)
% refinement for all sides
  sinterfin=cell(nv,1);
  sinterdifffin=cell(nv,1);
  npanfin=zeros(nv,1);
    
  for i=1:nv
    [sinterfin{i},sinterdifffin{i},npanfin(i)]=cornerrefine0(sinter{i}, ...
		sinterdiff{i},npan(i),nsub);
  end
%    
%    
%    
  function [sinter,sinterdiff,npanfin]=cornerrefine0(sinter0,sinterdiff0, ...
						     npan,nsub)
% refine each corner panel dyadically into nsub panels for one side
  npanfin=npan+2*nsub;

  sinter=zeros(npanfin+1,1);
  sinter(1:npan+1)=sinter0;
  
  sinterdiff=zeros(npanfin,1);
  sinterdiff(1:npan)=sinterdiff0;
  % refinement on the left corner  
  for k=1:nsub
    sinter(3:end)=sinter(2:end-1);
    sinter(2)=(sinter(1)+sinter(2))/2;   
    sinterdiff(3:end)=sinterdiff(2:end-1);
    sinterdiff(2)=sinterdiff(1)/2;
    sinterdiff(1)=sinterdiff(1)/2;   
  end
  % refinement on the right corner
  sinter(end-nsub:end)=1-flipud(sinter(1:nsub+1));
  sinterdiff(end-nsub-1:end)=flipud(sinterdiff(1:nsub+2));
    
  function M1=LogCinit(sinterdiff,T,W,npan,N,iper)
% *** excluding diagonal entries ***
% *** Corrections to Logarithmic potential log(|tau-z|) ***   
% *** block-tri-diagonal output ***
% *** iper=0,1 (0 is open arc, 1 is closed contour) ***
    Ng=length(T);
  M1=zeros(N);
  A=ones(Ng);
  for k=2:Ng
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
  TMP=-log(abs(T(:,ones(1,Ng))-T(:,ones(1,Ng))'));
  TMP(1:Ng+1:Ng*Ng)=0;
  TMP=TMP+WfrakLinit(A,0,1,T)./W(:,ones(1,Ng))';
  for k=1:npan
    myind=(k-1)*Ng+1:k*Ng;
    M1(myind,myind)=TMP;
  end
% *** superdiagonal blocks (targets to the left) ***
  for k=kstart:npan
    myinds=(k-1)*Ng+1:k*Ng;
    km1=mod(k-2,npan)+1;
    mi=(km1-1)*Ng+1:km1*Ng;       
    alpha=sinterdiff(km1)/sinterdiff(k);          
    [TMP,accept,na]=WfrakLinit(A,-1-alpha,alpha,T);
    mi=mi(accept);
    for nj=1:Ng
      M1(mi,myinds(nj))=-log(abs(T(nj)+1+(1-T(accept))*alpha));       
    end
    M1(mi,myinds)=M1(mi,myinds)+TMP./W(:,ones(1,na))';
  end
% *** subdiagonal blocks (targets to the right) ***
  for k=1:kend
    myinds=(k-1)*Ng+1:k*Ng;
    kp1=mod(k,npan)+1;
    mi=(kp1-1)*Ng+1:kp1*Ng;       
    alpha=sinterdiff(kp1)/sinterdiff(k);               
    [TMP,accept,na]=WfrakLinit(A,1+alpha,alpha,T);
    mi=mi(accept);
    for nj=1:Ng
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
  Ng=length(T);
  accept=1:Ng;
  T=trans+mscale*T;
  accept=accept(abs(T)<2);
  na=length(accept);
  p=zeros(1,Ng+1);
  Q=zeros(Ng,na);
  c=(1-(-1).^(1:Ng))./(1:Ng);
  for j=1:na
    jj=accept(j);
    p(1)=log(abs((1-T(jj))/(1+T(jj))));
    p111=log(abs(1-T(jj)^2));
    for k=1:Ng
      p(k+1)=T(jj)*p(k)+c(k);
    end
    Q(1:2:Ng-1,j)=(p111-p(2:2:Ng))./(1:2:Ng-1);
    Q(2:2:Ng,j)=(p(1)-p(3:2:Ng+1))./(2:2:Ng);
  end
  WfrakL=Q.'/A;
%
%
%
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
%
%
%
  function Pbc=Pbcinit(IP)
  Pbc=blkdiag(IP,rot90(IP,2));
%
%
%
  function Pbc=Poper(IP)
  Ng=length(IP)/2;
  Pbc=blkdiag(eye(Ng),IP,rot90(IP,2),eye(Ng));
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
  ns=5;
  theta=(0:ns-1)*2*pi/ns;
  zs=1.5*exp(1i*theta');
    
  f = zeros(size(ztarg));
  if nargout >1, fz = zeros(size(ztarg)); end
    
  for i=1:ns
    zi=ztarg-zs(i);
    r1=abs(zi);
    z=lambda*r1;
    
    f = f+q(i)*(log(r1)+besselk(0,z)); 
    if nargout > 1
      fz =fz+q(i)*(zi./r1.^2-lambda*zi.*besselk(1,z)./r1);
    end
  end

  function [wcmpL,wcmpC,wcmpH,wcmpS]=wLCHSinit(a,b,ztg,zsc,nz,wzp,awzp,io)
% *** zt is target vector in transformed plane ***
  Ng=length(zsc);
  cc=(b-a)/2;
  zt=(ztg-(b+a)/2)/cc;
  zsctr=(zsc-(b+a)/2)/cc;

  c1=(1-(-1).^(1:Ng))./(1:Ng);
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
%
%
%  
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
  
