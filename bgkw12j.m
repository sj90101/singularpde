  function bgkw12j
% *****************************************************************************
% * Solve the BGKW equation using the RCIP method
% *    
% * Authors: Johan Helsing and Shidong Jiang
% *
% * This function solves the so-called BGKW equation, see
% * Section 8 in https://arxiv.org/abs/2102.03504
% *
% * Last modified on 02/28/2021.
% *
% *****************************************************************************
  close all
  format long
  format compact
% *** Tables from Jiang and Luo JCP 2016 p. 416-434 ***************************
  table3=[0.003   -1.490909702131201e-3   1.242445655299172e-1;
          0.01    -4.900405009657547e-3   1.225330275292623e-1;
          0.03    -1.413798601526842e-2   1.180147037188893e-1;
          0.1     -4.155607782558620e-2   1.057028408172292e-1;
	  0.3     -9.344983511356682e-2   8.560111699820618e-2;
	  1.0     -1.694625753368226e-1   5.804708735555459e-2;
	  2.0     -2.083322536749375e-1   4.281659776113917e-2;
	  3.0     -2.266437497658084e-1   3.489298506190833e-2;
	  5.0     -2.446632678455994e-1   2.627042060967383e-2;
	  7.0     -2.536943539674479e-1   2.147460412330841e-2;
	  10.0    -2.611624603488405e-1   1.714449048590649e-2;
	  30.0    -2.743853873277227e-1   8.043009085700258e-3;
	  100.0   -2.796682147138912e-1   3.226757181742400e-3];
  table2=[0.003    0.4978915352789693    0.9939398014207545;
          0.01     0.4930697807742208    0.9800810020043015;
          0.03     0.4800058682766829    0.9425456003244014;
          0.1      0.4412246409722421    0.8352857656469417;
          0.3      0.3672125695500504    0.6635300770297120;
          1.0      0.2518613399894732    0.4442284697467991;
          2.0      0.1852462993740218    0.3274745769375018;
          3.0      0.1504282444992075    0.2672070059395808;
          5.0      0.1126351880294592    0.2016944318181617;
          7.0      0.09171689613521435   0.1652086347955840;
          10.0     0.07292211299328496   0.1321955790519311;
          30.0     0.03381357342231838   0.06243844733659697;
          100.0    0.01343072948081877   0.02520098284529616];
  kval=table3(:,1); % values of the Knudsen number
  pref=table3(:,2); % values of the stress P_{xy}
  qref=table3(:,3); % values of the half-channel mass flow rate Q=\int_0^0.5\rho(x)dx
  u1ref=table2(:,2); % velocity at x=1;
% ****************************************************
  T16=Tinit16;
  W16=Winit16;
  Ng=16;
  [~,~,U,~]=legeexps(Ng); % Legendre expansion conversion matrix
  A=ones(16); % The Vandermonde matrix
  for k=2:16
    A(:,k)=T16.*A(:,k-1);
  end  
%
% *** Pbc and PWbc prolongation matrices ***
  [IP,IPW]=IPinit(T16,W16);
  Pbc  = IP;
  PWbc = IPW;
  PbcL = blkdiag(Pbc,eye(16));
  starS=1:16;
  starL=1:32;
  circL=33:48;
  sidiloc=[0.5;0.5;1];
  LogCloc=LogCinit(sidiloc,T16,W16,3,48,0);  
%
  nsub=41;
  nplot=length(kval);
  qval=zeros(nplot,1);
  pval=zeros(nplot,1);
  u1comp=zeros(nplot,1);
  for k=1:nplot
    lambda=kval(k); % Knudsen number
    npan=4;
    hmax=1.5*lambda;
    h=1/npan;
    disp(['Knudsen number = ',num2str(lambda), ';  npan = ', num2str(npan)])
% *** Panels, discretization points, and weights ***
    sinter0=linspace(-0.5,0.5,npan+1)';
    sinterdiff=ones(npan,1)/npan;
    [z0,w,np]=zinit(sinter0,sinterdiff,T16,W16,npan);
    LogC=LogCinit(sinterdiff,T16,W16,npan,np,0);
    AbsC=AbsCinit(T16,W16);
    starind1=1:32; % left end point at 0
    starind2=np-31:np; % right end point at 1
%
% *** The K_coa^\circ matrix is set up ***
    t1=cputime;
    if hmax>h
      Kcirc=MAinit(LogC,AbsC,lambda,z0,w,sinterdiff,np);
    else
      isrcomp=0;
      Kcirc=MAadapinit(lambda,z0,w,sinter0,T16,W16,Ng,npan,hmax,isrcomp,A,U);
    end
    dt=cputime-t1;
    disp(['Coarse matrix construction time = ', num2str(dt), ' seconds'])
    Kcirc(starind1,starind1)=zeros(32);
    Kcirc(starind2,starind2)=zeros(32);
%
    fcoa=f0(z0,lambda); % the singular right-hand side
    fcirccoa=fcoa;
    fcirccoa(starind1)=0;
    fcirccoa(starind2)=0;
%
% *** Recursion for the R matrix ***
    R=speye(np);
    [R(starind1,starind1),rf,Rstor,rfstor,Kstor,fibstor]=Rcomp(LogCloc, ...
	AbsC,lambda,T16,W16,Pbc,PWbc,starS,starL,circL,nsub,npan,hmax,A,U);
% R matrix at the right end point by symmetry
    R(starind2,starind2)=rot90(R(starind1,starind1),2);
%
% *** Solving main linear system for vtilde ***
    rfstarcoa1=rf; % singular part at the left end point
    rfstarcoa2=-rot90(rf,2); % singular part at the right end point by symmetry
    rhs=fcirccoa-Kcirc(:,starind1)*rfstarcoa1-Kcirc(:,starind2)*rfstarcoa2;
    t1=cputime;
    [vtilde,it]=myGMRESR(Kcirc,R,rhs,np,400,eps);
    dt=cputime-t1;
    disp(['Time on GMRES = ', num2str(dt), ' seconds'])
    disp(['GMRES iter = ',num2str(it)])
    disp ' '   
    %cond(eye(np)+Kcirc*R)
% *** Reconstruction of the solution on the fine grid.
    [sinterfin0,sinterfindiff,npanfin]=panelinit0(nsub,npan);
    [zfin0,~,npfin]=zinit(sinterfin0,sinterfindiff,T16,W16,npanfin);
%    
    pts2=16*(2+nsub);
% *** Reconstruction of vfin on Gamma^{\star} ***
    vt=vtilde(starind1);
    vfinloc=zeros(pts2,1);
    for i=nsub:-1:1
      Kcirc=Kstor(:,:,i);
      MAT=eye(48)+Kcirc;
      tmp=PbcL*vt;
      vt=tmp-Kcirc*SchurBanaLs(MAT,Rstor(:,:,i),tmp,starL,circL);
      myindL=(i+1)*16+(1:16);
      vfinloc(myindL)=vt(circL);
      vt=vt(starL);
    end
    vfinloc(1:32)=Rstor(:,:,1)*vt;
    vfin=[vfinloc;vtilde(33:np-32);-flipud(vfinloc)];
%
% *** Reconstruction of gfin on Gamma^{\star} ***
    gt=zeros(32,1);   % smooth part
    gfinloc=zeros(pts2,1);
    for i=nsub:-1:1
      Kcirc=Kstor(:,:,i);
      MAT=eye(48)+Kcirc;    
      rf=fibstor(:,i);
      rf(starL)=rfstor(:,i);
      tmp=PbcL*gt;
      gt=tmp-Kcirc*(SchurBanaLs(MAT,Rstor(:,:,i),tmp,starL,circL) ...
		      +SchurBanaLsrf(MAT,Rstor(:,:,i),rf,starL,circL));
      myindL=(i+1)*16+(1:16);
      gfinloc(myindL)=fibstor(circL,i)+gt(circL);
      gt=gt(starL);    
    end
    gfinloc(1:32)=rfstor(:,1)+Rstor(:,:,1)*gt;
    gfin=[gfinloc;zeros(np-64,1);-flipud(gfinloc)];
%
% *** Post processing ***
    rhofin=vfin+gfin;
    rhohat=R*vtilde;
    rhohat(starind1)=rhohat(starind1)+rfstarcoa1;
    rhohat(starind2)=rhohat(starind2)+rfstarcoa2;
    if lambda<=0.3
      % add x
      rhofin=rhofin+zfin0;
      rhohat=rhohat+z0;        
    end
% evaluate the half-channel mass flow rate
% by rhohat. npan must be even for this to be accurate.    
    qval(k)=w(np/2+1:end)'*rhohat(np/2+1:end);
% evaluate the stress P_{xy}
    pval(k)=stresscomp0(lambda,rhohat,z0,w,sinter0,sinterdiff,npan,np, ...
			T16,W16,A,U,hmax);
    %u1comp(k)=ones(1,Ng)*U*rhofin(npfin-15:end);
    u1comp(k)=(ones(1,Ng)/A)*rhofin(npfin-15:end);
  end
%
  disp ' '
  rerr=abs((qref-qval)./qref);
  disp('   k            Q value (2016)          Q value (now)             Error')  
  fprintf('%7.3f &  %0.15e &    %0.15e &    %0.1e \\\\\n', [kval qref qval rerr]')
  disp ' '
  disp(['maximum relative error in Q value = ',num2str(max(rerr))])
  disp ' '
%  
  disp ' '
  rerr=abs((pref-pval)./pref);
  disp('   k            P value (2016)            P value (now)             Error')  
  fprintf('%7.3f &  %0.15e &    %0.15e  &   %0.1e \\\\\n', [kval pref pval rerr]')
  disp ' '
  disp(['maximum relative error in P value = ',num2str(max(rerr))])
  disp ' '
%  
  disp ' '
  rerr=abs((u1ref-u1comp)./u1ref);
  disp('   k            u(1) (2016)             u(1) (now)                Error')  
  fprintf('%7.3f &  %0.15e &    %0.15e  &   %0.1e \\\\\n', [kval u1ref u1comp rerr]')
  disp ' '
  disp(['maximum relative error in u(1) = ',num2str(max(rerr))])

  function pval=stresscomp0(lambda,rhohat,z,w,sinter,sinterdiff, ...
			    npan,np,T16,W16,A,U,hmax)
% evaluate the stress 
% P_{xy}=-\left( 2/k*\int_0.5^1 J_0((x-0.5)/k)rho(x)dx
%               +J_1(1/2/k) \right)/sqrt(pi)
  Ng=length(T16);
  targ=0.0;
  h=sinterdiff(npan/2+1);  
  if h > hmax
    nlevel=ceil(log2(h/hmax));
    h0=h/2^nlevel;
    t1=sinter(npan/2+1);
    t2=sinter(npan/2+2);
    [sinter0,sinterdiff0,npan0]=rightdivide(t1,t2,h0);
    [z0,w0,~]=zinit(sinter0,sinterdiff0,T16,W16,npan0);
    P0= Pinit0(t1,t2,z0,Ng,A);
    %P0= Pinit(t1,t2,z0,Ng,U);
    r0 = abram0(abs(z0-targ)/lambda).*w0;
    for i=1:1 % change to 2 has little effect on accuracy
      ind=(i-1)*Ng+(1:Ng);
      a=sinter0(i);
      b=sinter0(i+1);
      wcmpL=wLogCinit1(a,b,targ,A);
      zd=abs(z0(ind)-targ)/lambda;
      [fs,flog]=j0parts(zd);
      dsd  = sinterdiff0(i)/2/lambda;
      scale= sinterdiff0(i)/2;
      r0(ind)=(fs+log(dsd)*flog).*w0(ind)+scale*flog.*wcmpL;
    end
    pval=r0'*P0*rhohat(np/2+(1:Ng));
    ind=np/2+Ng+1:np;
    pval=pval+abram0(abs(z(ind)-targ)/lambda)'*(w(ind).*rhohat(ind));
    pval=-(pval*2/lambda+abram1(1/2/lambda))/sqrt(pi);
  else
    ising=1; % chenge to 2 has little effect on accuracy
    pval=0;
    for i=1:ising
      ind=(npan/2+i-1)*Ng+(1:Ng);
      a=sinter(npan/2+i);
      b=sinter(npan/2+i+1);
      wcmpL=wLogCinit1(a,b,targ,A);
      zd=abs(z(ind)-targ)/lambda;
      [fs,flog]=j0parts(zd);
      dsd=sinterdiff(npan/2+i)/2/lambda;
      scale=sinterdiff(npan/2+i)/2;
      rval=(fs+log(dsd)*flog).*w(ind)+scale*flog.*wcmpL;
      pval=pval+rval'*rhohat(ind);
    end
    ind=np/2+ising*Ng+1:np;
    pval=pval+abram0(abs(z(ind)-targ)/lambda)'*(w(ind).*rhohat(ind));
    pval=-(pval*2/lambda+abram1(1/2/lambda))/sqrt(pi);
  end
  
  function fout=f0(z0,lambda)
% the right hand side  
  if lambda >0.3
    fout= (abram0((0.5-z0)/lambda)-abram0((0.5+z0)/lambda))/(2*sqrt(pi));
  else
    fout=-(abram1((0.5-z0)/lambda)-abram1((0.5+z0)/lambda))*lambda/sqrt(pi);
  end
    
  function fout=f(z,lambda)
% the right hand side used in Rcomp, the left endpoint is shifted to 0.
  if lambda >0.3
    fout=(abram0((1-z)/lambda)-abram0(z/lambda))/(2*sqrt(pi));
  else
    fout=-(abram1((1-z)/lambda)-abram1(z/lambda))*lambda/sqrt(pi);
  end
  
  function [R,rf,Rstor,rfstor,Kstor,fibstor]=Rcomp(LogC,AbsC,lambda, ...
			T16,W16,Pbc,PWbc,starS,starL,circL,nsub,npan,hmax,A,U)
  Rstor=zeros(32,32,nsub);
  rfstor=zeros(32,nsub);
  Kstor=zeros(48,48,nsub);
  fibstor=zeros(48,nsub);  
  for level=1:nsub
    [z,w,sinterloc,sidi]=zlocinit(T16,W16,nsub,level,npan);
    floc=f(z,lambda);
    if hmax>sidi(3)
      K=MAinit(LogC,AbsC,lambda,z,w,sidi,48);
    else
      ifrcomp=1;
      K=MAadapinit(lambda,z,w,sinterloc,T16,W16,16,3,hmax,ifrcomp,A,U);
    end
    MAT=eye(48)+K;
    K(starL,starL)=zeros(32);
    if level==1
      R=inv(MAT(starL,starL));
      rf=R*floc(starL);
    end
    Rstor(:,:,level)=R;  
    rfstor(:,level)=rf;
    Kstor(:,:,level)=K;  
    fibstor(:,level)=floc;
    rf=SchurBanarf(PWbc,MAT,R,rf,floc,starS,starL,circL);
    R=SchurBana(Pbc,PWbc,MAT,R,starS,starL,circL);
  end

  function A=SchurBana(P,PW,K,A,starS,starL,circL)
  circS=17:32;
  VA=K(circL,starL)*A;
  PTA=PW'*A;
  PTAU=PTA*K(starL,circL);
  DVAUI=inv(K(circL,circL)-VA*K(starL,circL));
  DVAUIVAP=DVAUI*(VA*P);
  A(starS,starS)= PTA*P+PTAU*DVAUIVAP;
  A(circS,circS)= DVAUI;
  A(circS,starS)=-DVAUIVAP;
  A(starS,circS)=-PTAU*DVAUI;

  function A=SchurBanaRfstar(PW,K,A,starS,circS,starL,circL)
  AU=A*K(starL,circL);
  DVAUIV=(K(circL,circL)-K(circL,starL)*AU)\K(circL,starL);
  A(starS,:)= PW'+PW'*AU*DVAUIV;
  A(circS,:)=-DVAUIV;
  
  function rfout=SchurBanarf(PW,K,A,rf,f,starS,starL,circL)
  circS=17:32;
  AU=A*K(starL,circL);
  DVAU=K(circL,circL)-K(circL,starL)*AU;
  tmp=DVAU\(f(circL)-K(circL,starL)*rf);
  rfout=zeros(size(rf));
  rfout(starS)=PW'*(rf-AU*tmp);
  rfout(circS)=tmp;
  
  function x=SchurBanaLs(K,A,rhs,starL,circL)
  AU=A*K(starL,circL);
  DVAU=K(circL,circL)-K(circL,starL)*AU;
  Arhs=A*rhs(starL);
  x=zeros(size(rhs));
  tmp=DVAU\(rhs(circL)-K(circL,starL)*Arhs);
  x(starL)=Arhs-AU*tmp;
  x(circL)=tmp;
  
  function x=SchurBanaLsrf(K,A,rhs,starL,circL)
  AU=A*K(starL,circL);
  DVAU=K(circL,circL)-K(circL,starL)*AU;
  x=zeros(size(rhs));
  tmp=DVAU\(rhs(circL)-K(circL,starL)*rhs(starL));
  x(starL)=rhs(starL)-AU*tmp;
  x(circL)=tmp;
  
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
      y=triu(H(1:it,1:it))\s(1:it);             
      x=fliplr(V(:,1:it))*flipud(y);
      trueres=norm(x+A*(R*x)-b)/bnrm2;
      disp(['GMRES residual         = ',num2str(trueres)])
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

  function M1=MAinit(LogC,AbsC,lambda,z,w,sinterdiff,N)
% *** the kernel K(x,y)=J_{-1}(|x-y|/lambda) ***   
  W=w(:,ones(1,N))';
  zdiff=abs(z-z.')/lambda;
  M1=abramm1(zdiff).*W;  
% fix the logarithimic singularity  
  myind=find(LogC);
  ML=jm1log(zdiff(myind));
  M1(myind)=M1(myind)+ML.*LogC(myind).*W(myind);
  [MS0, ~, ML0]=jm1parts(0);  
  for m=1:N
    dsd=sinterdiff(fix((m-1)/16)+1)/2;
    M1(m,m)=( MS0 + ML0*(LogC(m,m)+log(dsd)-log(lambda)) )*w(m);
  end
% fix the absolute value singularity, only diagonal blocks 
  npan = N/16;
  for k=1:npan
    ind=(k-1)*16+1:k*16;
    MA=jm1abs(abs(z(ind)-z(ind).')/lambda);
    scale=sinterdiff(k)/2/lambda;
    M1(ind,ind)= M1(ind,ind)+scale*MA.*AbsC.*W(ind,ind);
  end    
  M1=-M1/(lambda*sqrt(pi));

  function M1=AbsCinit(T,W)
% *** Corrections to absolute singularity |x-y| ***   
  A=ones(16);
  for k=2:16
    A(:,k)=T.*A(:,k-1);
  end
  fint=zeros(16);
% fint(i,j)=\int_{-1}^1 |T_i-y| y^{j-1}dy
  for j=1:16
    fint(:,j)=2/(j*(j+1))*T.^(j+1)+(1-(-1)^j)/(j+1)-(1+(-1)^j)*T/j;
  end
  M1=fint/A;
  TMP=-abs(T-T');
  M1=TMP+M1./W(:,ones(1,16))';
  
  function M1=LogCinit(sinterdiff,T,W,npan,N,iper)
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
  TMP=-log(abs(T-T'));
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

  function [WfrakL,accept,na]=WfrakLinit(A,trans,mscale,T)
% *** T is target vector, sources on canonical interval ***
  accept=1:16;
  T=trans+mscale*T;
  accept=accept(abs(T)<2);
  na=length(accept);
  p=zeros(17,1);
  Q=zeros(16,na);
  c=(1-(-1).^(1:16))./(1:16);
  for j=1:na
    jj=accept(j);
    p(1)=log(abs((1-T(jj))/(1+T(jj))));
    p111=log(abs(1-T(jj)^2));
    for k=1:16
      p(k+1)=T(jj)*p(k)+c(k);
    end
    Q(1:2:15,j)=(p111-p(2:2:16))./(1:2:15)';
    Q(2:2:16,j)=(p(1)-p(3:2:17))./(2:2:16)';
  end
  WfrakL=Q.'/A;
  
  function M1=MAadapinit(lambda,z,w,sinter,T,W,Ng,npan,hmax,ifrcomp,A,U)
% *** the kernel K(x,y)=J_{-1}(|x-y|/lambda) ***
  N=Ng*npan;
  M1=zeros(N);
  W1=w(:,ones(1,N))';
  for i=1:npan
    im1=max(1,i-1);
    ip1=min(npan,i+1);
    ind1=(i-1)*Ng+(1:Ng);
    for j=setdiff(1:npan,im1:ip1)
      ind2=(j-1)*Ng+(1:Ng);
      zd=abs(z(ind1)-z(ind2).')/lambda;
      M1(ind1,ind2)=abramm1(zd).*W1(ind1,ind2);
    end
  end
  inds=1:Ng;
  if ifrcomp
    kstart=1;
    kend=npan;
  else
    kstart=2;
    kend=npan-1;
  end
  for k=kstart:kend
    kstart=(k-1)*Ng;
    km1=max(1,k-1);
    kp1=min(npan,k+1);
    t1=sinter(k);
    t2=sinter(k+1);
    for i=1:Ng
      targ=z(kstart+i);
      % self coarse panel
      [sinter0,sinterdiff0,npan0,m]=selfdivide(t1,t2,targ,hmax);
      [z0,w0,~]=zinit(sinter0,sinterdiff0,T,W,npan0);
      P0=Pinit(t1,t2,z0,Ng,U);
      %P0=Pinit0(t1,t2,z0,Ng,A);
      r0 = abramm1(abs(z0-targ)/lambda).*w0;
      % self fine panel
      ind0=(m-1)*Ng+(1:Ng);
      a=sinter0(m);
      b=sinter0(m+1);
      wcmpL=wLogCinit0(a,b,targ,A); % logarithmic singularity correction
      %wcmpL=wLogCinit(a,b,targ,U); % logarithmic singularity correction
      wcmpA=wAbsCinit0(a,b,targ,A); % absolute value singularity correction
      %wcmpA=wAbsCinit(a,b,targ,U); % absolute value singularity correction
      [fs,fabs,flog]=jm1parts(abs(z0(ind0)-targ)/lambda);
      dsd    =sinterdiff0(m)/2/lambda;
      scale =(sinterdiff0(m)/2)^2/lambda;
      lscale =sinterdiff0(m)/2;
      r0(ind0)=(fs + log(dsd)*flog).*w0(ind0) ...
		+ lscale*flog.*wcmpL + scale*fabs.*wcmpA;
      mlr=[];
      if m>1 && m<npan0
        mlr=[m-1 m+1];
      elseif m==1 && m<npan0
        mlr=m+1;
      elseif m==npan0 && m>1
        mlr=m-1;
      end
      % left and/or right fine panels on the self coarse panel
      for j=mlr
    	indlr=(j-1)*Ng+(1:Ng);
        a=sinter0(j);
        b=sinter0(j+1);
        cc=(b-a)/2;
        zt=(targ-(b+a)/2)/cc;
        if abs(zt)<2
          %wcmpL=wLogCinit(a,b,targ,U); % logarithmic singularity correction
          wcmpL=wLogCinit0(a,b,targ,A); % logarithmic singularity correction
          zd=abs(z0(indlr)-targ)/lambda;
          [fs,fabs,flog]=jm1parts(zd);
          dsd  =sinterdiff0(j)/2/lambda;
          scale=sinterdiff0(j)/2;
          r0(indlr)=(fs + fabs.*zd+log(dsd)*flog).*w0(indlr) ...
		    + scale*flog.*wcmpL;
    	end
      end
      M1(kstart+i,kstart+inds)=r0'*P0;
      if km1 < k % left corase panel
        a=sinter(km1);
        b=sinter(km1+1);
        h0=sinterdiff0(1);
        [sinterm1,sinterdiffm1,npanm1]=leftdivide(a,b,h0);
        [zm1,wm1,~]=zinit(sinterm1,sinterdiffm1,T,W,npanm1);
        Pm1=Pinit(a,b,zm1,Ng,U);
        %Pm1=Pinit0(a,b,zm1,Ng,A);
        rm1 = abramm1(abs(zm1-targ)/lambda).*wm1;
        if m==1 % left fine panel on the left coarse panel
          indl=(npanm1-1)*Ng+(1:Ng);
          a=sinterm1(npanm1);
          b=sinterm1(npanm1+1);
          cc=(b-a)/2;
	      zt=(targ-(b+a)/2)/cc;
          if (abs(zt)<2) % corrections for the logarithmic singularity
            %wcmpL=wLogCinit(a,b,targ,U);
            wcmpL=wLogCinit0(a,b,targ,A); % logarithmi singularity correction
            zd=abs(zm1(indl)-targ)/lambda;
            [fs,fabs,flog]=jm1parts(zd);
            dsd  =sinterdiffm1(npanm1)/2/lambda;
            scale=sinterdiffm1(npanm1)/2;
            rm1(indl)=(fs + fabs.*zd+log(dsd)*flog).*wm1(indl) ...
                      + scale*flog.*wcmpL;
          end
        end
        M1(kstart+i,(km1-1)*Ng+inds)=rm1'*Pm1;
      end
      if kp1 > k % right coarse panel
        a=sinter(kp1);
        b=sinter(kp1+1);
        h0=sinterdiff0(end);
        [sinterp1,sinterdiffp1,npanp1]=rightdivide(a,b,h0);
        [zp1,wp1,~]=zinit(sinterp1,sinterdiffp1,T,W,npanp1);
        Pp1= Pinit(a,b,zp1,Ng,U);
        %Pp1=Pinit0(a,b,zp1,Ng,A);
        rp1 = abramm1(abs(zp1-targ)/lambda).*wp1;
        if m==npan0 % right fine panel on the right coarse panel
	  indr=1:Ng;
          a=sinterp1(1);
          b=sinterp1(2);
	  cc=(b-a)/2;
	  zt=(targ-(b+a)/2)/cc;
          if (abs(zt)<2) % corrections for the logarithmic singularity
            %wcmpL=wLogCinit(a,b,targ,U);
            wcmpL=wLogCinit0(a,b,targ,A); % logarithmi singularity correction
            zd=abs(zp1(indr)-targ)/lambda;
            [fs,fabs,flog]=jm1parts(zd);
            dsd  =sinterdiffp1(1)/2/lambda;
            scale=sinterdiffp1(1)/2;
            rp1(indr)=(fs + fabs.*zd+log(dsd)*flog).*wp1(indr) ...
                      + scale*flog.*wcmpL;
          end
        end
        M1(kstart+i,(kp1-1)*Ng+inds)=rp1'*Pp1;
      end
    end
  end
  M1=-M1/(lambda*sqrt(pi));

  function Ploc=Pinit(a,b,z,Ng,U)
  cc=(b-a)/2;
  zt=(z-(b+a)/2)/cc;
  Ploc=legepols(zt,Ng-1)*U;
      
  function Ploc=Pinit0(a,b,z,Ng,A)
  AA=ones(length(z),Ng);
  cc=(b-a)/2;
  zt=(z-(b+a)/2)/cc;
  for k=2:16
    AA(:,k)=AA(:,k-1).*zt;   
  end
  Ploc=AA/A;
  
  function wcmpA=wAbsCinit0(a,b,ztg,A)
% *** t is the target point on the standard interval [-1,1] ***
  cc=(b-a)/2;
  t=(ztg-(b+a)/2)/cc;
  f=zeros(1,16);
  for j=1:16  % f(k+1) = \int_{-1}^1 |t-x|x^k dx
    f(j)=2/(j*(j+1))*t.^(j+1)+(1-(-1)^j)/(j+1)-(1+(-1)^j)*t/j;
  end
  wcmpA=(f/A)';
  
  function wcmpA=wAbsCinit(a,b,ztg,U)
% *** t is the target point on the standard interval [-1,1] ***
  Ng=16;
  cc=(b-a)/2;
  t=(ztg-(b+a)/2)/cc;
  p=legepols(t,Ng+1);
  f=zeros(1,Ng); % *** f(k+1) = \int_{-1}^1 |t-x|p_k(x)dx
  for j=1:2
    f(j)=2/(j*(j+1))*t.^(j+1)+(1-(-1)^j)/(j+1)-(1+(-1)^j)*t/j;
  end
  for j=3:Ng
    f(j)=2/(2*j-1)*((p(j+2)-p(j))/(2*j+1)+(p(j-2)-p(j))/(2*j-3));
  end
  wcmpA=(f*U)';
  
  function wcmpL=wLogCinit0(a,b,ztg,A)
% *** zt is target point on the standard interval [-1,1] ***
  cc=(b-a)/2;
  zt=(ztg-(b+a)/2)/cc;
  p=zeros(17,1);
  Q=zeros(16,1); % Q(k+1)=\int_{-1}^1 z^k*log(z-zt) dz
  c=(1-(-1).^(1:16))./(1:16);
  p(1)=log(abs((1-zt)/(1+zt)));
  p111=log(abs(1-zt.^2));
  for k=1:16
    p(k+1)=zt*p(k)+c(k);
  end
  Q(1:2:15)=(p111-p(2:2:16))./(1:2:15)';
  Q(2:2:16)=(p(1)-p(3:2:17))./(2:2:16)';
  wcmpL=A'\Q;
  
  function wcmpL=wLogCinit1(a,b,ztg,A)
% *** zt is target point on the standard interval [-1,1] ***
  cc=(b-a)/2;
  zt=(ztg-(b+a)/2)/cc;
  if abs(zt+1)>eps*100
    wcmpL=wLogCinit0(a,b,ztg,A);
  else
    Q=zeros(16,1); % Q(k+1)=\int_{-1}^1 z^k*log(z+1) dz
    c=(1-(-1).^(1:16))./(1:16);
    Q(1)=2*log(2)-2;
    for k=1:15
      Q(k+1)=(2*log(2)-c(k+1)-k*Q(k))/(k+1);
    end
    wcmpL=A'\Q;
  end
  
  function wcmpL=wLogCinit(a,b,ztg,U)
% *** zt is target point on the standard interval [-1,1] ***
  Ng=16;
  cc=(b-a)/2;
  zt=(ztg-(b+a)/2)/cc;
  f=zeros(Ng+1,1); % f(k+1)=\int_{-1}^1 p_k(z)/(z-zt) dz 
  q=zeros(1,Ng); % q(k+1)=\int_{-1}^1 p_k(z)*log(z-zt) dz  
  f(1)=log(abs((1-zt)/(1+zt)));
  f(2)=2+zt*f(1);
  for k=1:Ng-1
    f(k+2)  =( (2*k+1)*zt*f(k+1)-k*f(k) )/(k+1);
  end
  q(1)=2*log(abs(1-zt))-(zt+1)*f(1)-2;
  for k=1:Ng-1
    q(k+1)=(f(k)-f(k+2))/(2*k+1);
  end
  wcmpL=(q*U)';
  
  function [sinter,sinterdiff,npanfin]=panelinit0(nsub,npan)
  npanfin=npan+2*nsub;
  sinter=zeros(npanfin+1,1);
  sinter(1:npan+1)=linspace(-0.5,0.5,npan+1);
  sinterdiff=ones(npan+2*nsub,1)/npan;
  for k=1:nsub
    sinter(3:end)=sinter(2:end-1);
    sinter(2)=(sinter(1)+sinter(2))/2;   
    sinterdiff(3:end)=sinterdiff(2:end-1);
    sinterdiff(2)=sinterdiff(1)/2;
    sinterdiff(1)=sinterdiff(1)/2;   
  end
  sinter(end-nsub:end)=-flipud(sinter(1:nsub+1));
  sinterdiff(end-nsub-1:end)=flipud(sinterdiff(1:nsub+2));
  
  function [z,w,np,s]=zinit(sinter,sinterdiff,T,W,npan)
  np=16*npan;
  s=zeros(np,1);
  w=zeros(np,1);
  for k=1:npan
    myind=(k-1)*16+1:k*16;
    sdif=sinterdiff(k)/2;
    s(myind)=(sinter(k)+sinter(k+1))/2+sdif*T;
    w(myind)=W*sdif;
  end
  z=zfunc(s) ;

  function [z,w,sinterloc,sidi]=zlocinit(T,W,nsub,level,npan)
  denom=2^(nsub-level)*npan;
  s=[T/4+0.25;T/4+0.75;T/2+1.5]/denom;
  w=[W/4;W/4;W/2]/denom;
  sidi=[0.5;0.5;1]/denom;  
  sinterloc=[0; cumsum(sidi)];
  z=zfunc(s);

  function zout=zfunc(s)
  zout=s;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  functions for evaluating the Abramowitz functions
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function v = abram1(x)
  ab1f = [1.47285192577978807369d0; 0.10903497570168956257d0;
          -0.12430675360056569753d0; 0.306197946853493315d-2;
          -0.2218410323076511d-4; 0.6989978834451d-7;
	  -0.11597076444d-9; 0.11389776d-12; -0.7173d-16; 0.3d-19];
  ab1g = [0.39791277949054503528d0; -0.29045285226454720849d0;
          0.1048784695465363504d-1; -0.10249869522691336d-3;
          0.41150279399110d-6; -0.83652638940d-9;
          0.97862595d-12; -0.71868d-15; 0.35d-18];
  ab1h = [0.84150292152274947030d0; -0.7790050698774143395d-1;
          0.133992455878390993d-2; -0.808503907152788d-5;
          0.2261858281728d-7; -0.3441395838d-10;
          0.3159858d-13; -0.1884d-16; 0.1d-19];
  ab1as = [2.13013643429065549448d0; 0.6371526795218539933d-1;
	   -0.129334917477510647d-2; 0.5678328753228265d-4;
	   -0.279434939177646d-5; 0.5600214736787d-7;
	   0.2392009242798d-7; -0.750984865009d-8;
	   0.173015330776d-8; -0.36648877955d-9; 0.7520758307d-10;
	   -0.1517990208d-10; 0.301713710d-11; -0.58596718d-12;
	   0.10914455d-12; -0.1870536d-13; 0.262542d-14; -0.14627d-15;
	   -0.9500d-16; 0.5873d-16; -0.2420d-16; 0.868d-17; -0.290d-17;
	   0.93d-18; -0.29d-18; 0.9d-19; -0.3d-19; 0.1d-19];
  onerpi = 0.56418958354775628695d0;
  rt3bpi = 0.97720502380583984317D+00;
  xlow = 1.11023d-16;
  xlow1 = 1.490116d-8;
  lnxmin = -708.3964d0;
  ntermf=9;
  nterma=27;
  ab1as(1) = ab1as(1)/2;
  ab1f(1) = ab1f(1)/2;
  ab1g(1) = ab1g(1)/2;
  ab1h(1) = ab1h(1)/2;
  v=zeros(size(x));
  if sum(x<0)>0
    error('Abram1 - Fatal error: argument x<0');
  end
  v(x<xlow) = 0.5;
  ind = find(x>=xlow & x<xlow1);
  if ~isempty(ind)
    v(ind)=(1 - x(ind)/onerpi - x(ind).* x(ind).* log(x(ind)))*0.5;
  end
  ind = find(x>=xlow1 & x<=2);
  if ~isempty(ind)
    t=x(ind).*x(ind)/2-1;
    pols=chebpols(t,ntermf);
    fval = pols*ab1f;
    gval = pols(:,1:end-1)*ab1g;
    hval = pols(:,1:end-1)*ab1h;
    v(ind) = fval-x(ind).*(gval/onerpi+x(ind).*log(x(ind)).*hval);
  end
  ind = find(x>2);
  if ~isempty(ind)
    nu = 3*(x(ind)/2).^(2/3);
    t = 6./nu -1;
    pols = chebpols(t,nterma);
    asval = pols*ab1as;
    asln = log(asval.*sqrt(nu/3)/rt3bpi) - nu;
    ind2 = find(asln>=lnxmin);
    u=zeros(size(asln));
    u(ind2)=exp(asln(ind2));
    v(ind)=u;
  end
	  
  function v = abram0(x)
% evaluate the Abramowitz function of order zero J_0(x), 
% translation of the Fortran code by McLeod 
% separates the logarithmic part when x is small
  ab0f=[-0.68121927093549469816d0; -0.78867919816149252495d0;
	0.5121581776818819543d-1; -0.71092352894541296d-3;
	0.368681808504287d-5; -0.917832337237d-8; 0.1270202563d-10;
	-0.1076888d-13; 0.599d-17];
  ab0g=[-0.60506039430868273190d0; -0.41950398163201779803d0;
	0.1703265125190370333d-1; -0.16938917842491397d-3;
	0.67638089519710d-6; -0.135723636255d-8;
	0.156297065d-11; -0.112887d-14; 0.55d-18];
  ab0h=[1.38202655230574989705d0; -0.30097929073974904355d0;
	0.794288809364887241d-2; -0.6431910276847563d-4;
	0.22549830684374d-6; -0.41220966195d-9;
	0.44185282d-12; -0.30123d-15; 0.14d-18];
  ab0as=[1.97755499723693067407d+0; -0.1046024792004819485d-1;
	 0.69680790253625366d-3; -0.5898298299996599d-4; 0.577164455305320d-5;
	 -0.61523013365756d-6; 0.6785396884767d-7; -0.723062537907d-8;
	 0.63306627365d-9; -0.989453793d-11; -0.1681980530d-10; 0.673799551d-11;
	 -0.200997939d-11; 0.54055903d-12; -0.13816679d-12; 0.3422205d-13;
	 -0.826686d-14; 0.194566d-14; -0.44268d-15; 0.9562d-16; -0.1883d-16;
	 0.301d-17; -0.19d-18; -0.14d-18; 0.11d-18; -0.4d-19; 0.2d-19; -0.1d-19];
  gval0 = 0.13417650264770070909D+00;
  lnxmin = -708.3964D+00;
  nterma = 27;
  ab0as(1) = ab0as(1)/2;
  ab0f(1) = ab0f(1)/2;
  ab0g(1) = ab0g(1)/2;
  ab0h(1) = ab0h(1)/2;
  ntermf = 8;
  onerpi = 0.56418958354775628695D+00;
  rt3bpi = 0.97720502380583984317D+00;
  rtpib2 = 0.88622692545275801365D+00;
  xlow1 = 1.490116D-08;
  v=zeros(size(x));
  if sum(x<0)>0
    error('Abram0 - Fatal error: argument x<0');
  end
  v(x==0)=rtpib2;
  ind=find(x~=0 & x<xlow1);
  if ~isempty(ind)
    v(ind)=rtpib2+x(ind).*(log(x(ind))-gval0);
  end
  ind=find(x>=xlow1 & x<=2);
  if ~isempty(ind)
    t=x(ind).*x(ind)/2-1;
    pols=chebpols(t,ntermf);
    fval=pols*ab0f;
    gval=pols*ab0g;
    hval=pols*ab0h;
    v(ind)=fval/onerpi+x(ind).*(log(x(ind)).*hval-gval);
  end
  ind=find(x>2);
  if ~isempty(ind)
    nu=3*(x(ind)/2).^(2/3);
    t=6./nu-1;
    pols=chebpols(t,nterma);
    asval=pols*ab0as;
    asln=log(asval/rt3bpi)-nu;
    ind2=find(asln>=lnxmin);
    u=zeros(size(asln));
    u(ind2)=exp(asln(ind2));
    v(ind)=u;
  end

  function [fsmooth,flog]=j0parts(x)
% kernel split for the Abramowitz function of order zero J_0(x), 
  ab0f = [-0.68121927093549469816d0; -0.78867919816149252495d0;
	  0.5121581776818819543d-1; -0.71092352894541296d-3;
	  0.368681808504287d-5; -0.917832337237d-8; 0.1270202563d-10;
	  -0.1076888d-13; 0.599d-17];
  ab0g = [-0.60506039430868273190d0; -0.41950398163201779803d0;
	  0.1703265125190370333d-1; -0.16938917842491397d-3;
	  0.67638089519710d-6; -0.135723636255d-8;
	  0.156297065d-11; -0.112887d-14; 0.55d-18];
  ab0h = [1.38202655230574989705d0; -0.30097929073974904355d0;
	  0.794288809364887241d-2; -0.6431910276847563d-4;
	  0.22549830684374d-6; -0.41220966195d-9;
	  0.44185282d-12; -0.30123d-15; 0.14d-18];
  a = [-2.0000000000000000000000000000000000000000000000000D0;
       +3.3333333333333333333333333333333333333333333333333D-1;
       -8.3333333333333333333333333333333333333333333333333D-3;
       +6.6137566137566137566137566137566137566137566137566D-5;
       -2.2964432686654908877131099353321575543797766019988D-7;
       +4.1753513975736197958420180642402864625086847309070D-10;
       -4.4608455102282262776089936583763744257571418065245D-13;
       +3.0345887824681811412306079308682819222837699364112D-16;
       -1.3945720507666273626978896741122619128142325075419D-19;
       +4.5307733943035326923258274012744051748350633773292D-23;
       -1.0787555700722696886490065241129536130559674707927D-26];
  b = [-1.7724538509055160272981674833411451827975494561224D0;
       +2.6835300529540141818046372975279270687352199218023D-1;
       +1.7724538509055160272981674833411451827975494561224D0;
       -4.8916994532701134747452173273657656225669810980782D-1;
       -9.8469658383639779294342637963396954599863858673466D-2;
       +1.8062581966508617020196376651747747389750786078529D-2;
       +1.3129287784485303905912351728452927279981847823129D-3;
       -1.7484790424414750217590976682541044161681777814907D-4;
       -6.6986162165741346458736488410474118775417590934331D-6;
       +6.9003789666065488683377121580855861174433144031200D-7;
       +1.6539793127343542335490490965549165129732738502304D-8;
       -1.3760791254942414610313518269752785366263098109209D-9;
       -2.2782084197442895778912521991114552520293028240088D-11;
       +1.5788316272729113463256804298998502922367229103736D-12;
       +1.9257890276790275383696130170003848284271367912162D-14;
       -1.1376169807979124663099205130412199832946006968810D-15;
       -1.0698827931550152990942294538891026824595204395646D-17;
       +5.4843743696166248199349085478203605689667738841116D-19;
       +4.1133517614571906924038041287547200402134580529203D-21;
       -1.8559861431739840375659942328929793988883889143369D-22];
  ab0f(1) = ab0f(1)/2;
  ab0g(1) = ab0g(1)/2;
  ab0h(1) = ab0h(1)/2;
  ntermf = 8;
  onerpi = 0.56418958354775628695D+00;
  xlow1 = 1;
  fsmooth=zeros(size(x));
  flog=zeros(size(x));
  if sum(x<0 | x>2)>0
    error('j0parts - Fatal error: argument x<0 or x>2');
  end
  ind=find(x<xlow1);
  if ~isempty(ind)
    m=19;
    c1=b(m);
    for k=m-1:-1:1
      c1=c1.*x(ind)+b(k);
    end
    m2=9;
    x2=x(ind).^2;
    c2=a(m2);
    for k=m2-1:-1:1
      c2=c2.*x2+a(k);
    end
    fsmooth(ind) = -c1/2;
    flog(ind) = -c2.*x(ind)/2;
  end
  ind=find(x>=xlow1 & x<=2);
  if ~isempty(ind)
    t=x(ind).*x(ind)/2-1;
    pols=chebpols(t,ntermf);
    fval=pols*ab0f;
    gval=pols*ab0g;
    hval=pols*ab0h;
    fsmooth(ind)=fval/onerpi-x(ind).*gval;
    flog(ind)=x(ind).*hval;
  end
  
  function v = abramm1(x)
% evaluate the Abramowitz function of order -1, i.e., J_{-1}(x)
% separates the logarithmic part when x is small  
  ab0as=[1.97755499723693067407d+0; -0.1046024792004819485d-1;
         0.69680790253625366d-3; -0.5898298299996599d-4; 0.577164455305320d-5;
        -0.61523013365756d-6; 0.6785396884767d-7; -0.723062537907d-8;
         0.63306627365d-9; -0.989453793d-11; -0.1681980530d-10;
         0.673799551d-11; -0.200997939d-11; 0.54055903d-12;
        -0.13816679d-12; 0.3422205d-13; -0.826686d-14; 0.194566d-14;
        -0.44268d-15; 0.9562d-16; -0.1883d-16; 0.301d-17; -0.19d-18;
        -0.14d-18; 0.11d-18; -0.4d-19; 0.2d-19; -0.1d-19];
  ab2as=[2.46492325304334856893d0; 0.23142797422248905432d0;
	-0.94068173010085773d-3; 0.8290270038089733d-4; -0.883894704245866d-5;
         0.106638543567985d-5; -0.13991128538529d-6; 0.1939793208445d-7;
        -0.277049938375d-8; 0.39590687186d-9; -0.5408354342d-10; 0.635546076d-11;
        -0.38461613d-12; -0.11696067d-12; 0.6896671d-13; -0.2503113d-13;
         0.785586d-14; -0.230334d-14; 0.64914d-15; -0.17797d-15; 0.4766d-16;
        -0.1246d-16; 0.316d-17; -0.77d-18; 0.18d-18; -0.4d-19; 0.1d-19];
  ab0as(1) = ab0as(1)/2;
  ab2as(1) = ab2as(1)/2;
  gval0 = 0.86582349735229929091d0;
  lnxmin = -708.3964D+00;
  nterma0 = 27;
  rt3bpi = 0.97720502380583984317D+00;
  xlow = 2.22045d-16;
  xlow1 = 1.490116d-4;
  a=[-2.0000000000000000000000000000000000000000000000000D0;
     +1.0000000000000000000000000000000000000000000000000D0;
     -4.1666666666666666666666666666666666666666666666667D-2;
     +4.6296296296296296296296296296296296296296296296296D-4;
     -2.0667989417989417989417989417989417989417989417989D-6;
     +4.5928865373309817754262198706643151087595532039976D-9;
     -5.7990991632966941608916917558892867534842843484819D-12;
     +4.5518831737022717118459118963024228834256549046169D-15;
     -2.3707724863032665165864124459908452517841952628213D-18];
  b=[-1.7316469947045985818195362702472072931264780078198D0;
     +3.5449077018110320545963349666822903655950989122448D0;
     -1.1341765026477007090902318648763963534367609960901D0;
     -3.9387863353455911717737055185358781839945543469386D-1;
     +8.1979576499209751767648549925405403615420597059310D-2;
     +7.8775726706911823435474110370717563679891086938773D-3;
     -1.1577977635714663776652308016403069537515868809059D-3;
     -5.3588929732593077166989190728379295020334072747464D-5;
     +5.9806967430793448927326299487438117502610053026081D-6;
     +1.6539793127343542335490490965549165129732738502304D-7;
     -1.4719335240679294091760668290304035256638539447040D-8;
     -2.7338501036931474934695026389337463024351633888106D-10;
     +2.0078726603525024874472946222860416356501683654204D-11;
     +2.6961046387506385537174582238005387597979915077027D-13;
     -1.6760795833721868880525746902531471557190633459574D-14;
     -1.7118124690480244785507671262225642919352327033033D-16;
     +9.1839792232715994576195555638833867759620923522355D-18];
  a = a/2;
  b = b/2;
  v=zeros(size(x));
  if sum(x<0)>0
    error('Abram - Fatal error: argument x<0');
  end
  v(x==0)=0; % fake value, should blow up when x is equal to 0
  ind=find(x~=0 & x<xlow);
  if ~isempty(ind)
    v(ind)=-log(x(ind))-gval0;
  end
  ind=find(x>=xlow & x<xlow1);
  if ~isempty(ind)
    x2=x(ind).^2;
    x4=x2.^2;
    x6=x2.*x4;
    x8=x4.*x4;    
    feven=b(1)+b(3)*x2+b(5)*x4+b(7)*x6+b(9)*x8;
    fodd=b(2)+b(4)*x2+b(6)*x4+b(8)*x6;
    v(ind)=log(x(ind)).*(-1+x2/2-x4/48+a(4)*x6+a(5)*x8)+feven+x(ind).*fodd;
  end
  ind=find(x>=xlow1 & x<=2);
  if ~isempty(ind)
    [feven,fodd,flog] =jm1parts(x(ind));
    v(ind) = flog.*log(x(ind))+feven+x(ind).*fodd;
  end
  ind=find(x>2);
  if ~isempty(ind)
    nu=3*(x(ind)/2).^(2/3);
    t=6./nu-1;
    pols=chebpols(t,nterma0);
    as0val=pols*ab0as; 
    as2val=pols(:,1:end-1)*ab2as; 
    asln=log((as2val.*nu*2/3-as0val)/rt3bpi./x(ind))-nu;
    ind2=find(asln>=lnxmin);
    u=zeros(size(asln));
    u(ind2)=exp(asln(ind2));
    v(ind)=u;
  end

  function [feven,fodd,flog] =jm1parts(x)
% use quadruple precision fortran code to generate
% Chebyshev expansion coefficients of 
% odd, even, and log parts of J_{-1}(x) for x \in [0,10]
%
% That is, J_{-1}(x) = feven + fodd * |x| + flog * log(x)
%
% Evaluation of J_{-1}(x) by this kernel splitting gradually
% loses precision as x increases due to cancellation, even
% though in quadruple precision it is accurate to double precision.
%
% Using double precision arithmetic, the precision is only
% 1e-11 for x>7.
%
  ceven0 = [0.9339121500154000E+01;  0.9467951631766124E+01;
	    -.6909142551630843E+01;  -.4887443926742919E+01;
	    0.1187511775377600E+01;  -.9348289247701976E-01;
	    0.3627591153380088E-02;  -.8180961654295485E-04;
	    0.1178872899827834E-05;  -.1156887976958917E-07;
	    0.8100517403852717E-10;  -.4194838224358275E-12;
	    0.1653254193882329E-14;  -.5076725423076273E-17;
	    0.1238733627437446E-19];
  codd0 = [0.2762625360277227E+00;  0.3936070023845599E+00;
           0.1418765951993083E+01;  -.4271154385223360E+00;
	   0.4189464696930939E-01;  -.1968910398647824E-02;
	   0.5248267880174669E-04;  -.8767988412756958E-06;
	   0.9823332711469947E-08;  -.7755067024809453E-10;
	   0.4481315148667998E-12;  -.1953791457550506E-14;
	   0.6588109029923584E-17;  -.1753990524488704E-19];
  clog0 = [-.5001323317427930E+01;  -.7493821150778837E+01;
	   -.5102067942811160E+00;  0.2508222355964408E+01;
	   -.4427800206422171E+00;  0.3019820668279420E-01;
	   -.1067673263388026E-02;  0.2247303511677814E-04;
	   -.3065868576163888E-06;  0.2875538347259317E-08;
	   -.1937394988871976E-10;  0.9702899915931817E-13;
	   -.3712968268530488E-15;  0.1110500759672965E-17;
	   -.2625105260334293E-20];
  a=[-2.0000000000000000000000000000000000000000000000000D0;
     +1.0000000000000000000000000000000000000000000000000D0;
     -4.1666666666666666666666666666666666666666666666667D-2;
     +4.6296296296296296296296296296296296296296296296296D-4;
     -2.0667989417989417989417989417989417989417989417989D-6;
     +4.5928865373309817754262198706643151087595532039976D-9;
     -5.7990991632966941608916917558892867534842843484819D-12;
     +4.5518831737022717118459118963024228834256549046169D-15;
     -2.3707724863032665165864124459908452517841952628213D-18];
  b=[-1.7316469947045985818195362702472072931264780078198D0;
     +3.5449077018110320545963349666822903655950989122448D0;
     -1.1341765026477007090902318648763963534367609960901D0;
     -3.9387863353455911717737055185358781839945543469386D-1;
     +8.1979576499209751767648549925405403615420597059310D-2;
     +7.8775726706911823435474110370717563679891086938773D-3;
     -1.1577977635714663776652308016403069537515868809059D-3;
     -5.3588929732593077166989190728379295020334072747464D-5;
     +5.9806967430793448927326299487438117502610053026081D-6;
     +1.6539793127343542335490490965549165129732738502304D-7;
     -1.4719335240679294091760668290304035256638539447040D-8;
     -2.7338501036931474934695026389337463024351633888106D-10;
     +2.0078726603525024874472946222860416356501683654204D-11;
     +2.6961046387506385537174582238005387597979915077027D-13;
     -1.6760795833721868880525746902531471557190633459574D-14;
     -1.7118124690480244785507671262225642919352327033033D-16;
     +9.1839792232715994576195555638833867759620923522355D-18];
  feven=zeros(size(x));
  fodd=zeros(size(x));
  flog=zeros(size(x));
  ind=find(x<4/5);
  if ~isempty(ind)
    a = a/2;
    b = b/2;
    x2 = x(ind).^2;
    m2 = 9;
    tmp = a(m2)*ones(size(x2));
    for i=m2-1:-1:1
      tmp=tmp.*x2+a(i);
    end
    flog(ind)=tmp;
    m = 17;
    tmp = b(m-1)*ones(size(x2));
    for i=m-3:-2:1
      tmp = tmp.*x2+b(i);
    end
    fodd(ind)=tmp;
    tmp = b(m)*ones(size(x2));
    for i=m-2:-2:1
      tmp = tmp.*x2+b(i);
    end
    feven(ind)=tmp;
  end
  ind=find(x>=4/5 & x<=10);
  if ~isempty(ind)
    neven = 14; 
    nodd  = 13; 
    L = 10/sqrt(2);
    t = x(ind).^2/L^2-1;
    pols=chebpols(t,neven);
    feven(ind)=pols*ceven0;
    fodd(ind)=pols(:,1:nodd+1)*codd0;
    flog(ind)=pols*clog0;
  end
  
  function fodd =jm1abs(x)
  codd0 = [0.2762625360277227E+00;  0.3936070023845599E+00; 
           0.1418765951993083E+01;  -.4271154385223360E+00; 
	   0.4189464696930939E-01;  -.1968910398647824E-02; 
	   0.5248267880174669E-04;  -.8767988412756958E-06;
	   0.9823332711469947E-08;  -.7755067024809453E-10;
	   0.4481315148667998E-12;  -.1953791457550506E-14;
	   0.6588109029923584E-17;  -.1753990524488704E-19];
  b=[-1.7316469947045985818195362702472072931264780078198D0; 
     +3.5449077018110320545963349666822903655950989122448D0; 
     -1.1341765026477007090902318648763963534367609960901D0; 
     -3.9387863353455911717737055185358781839945543469386D-1;
     +8.1979576499209751767648549925405403615420597059310D-2;
     +7.8775726706911823435474110370717563679891086938773D-3;
     -1.1577977635714663776652308016403069537515868809059D-3;
     -5.3588929732593077166989190728379295020334072747464D-5;
     +5.9806967430793448927326299487438117502610053026081D-6;
     +1.6539793127343542335490490965549165129732738502304D-7;
     -1.4719335240679294091760668290304035256638539447040D-8;
     -2.7338501036931474934695026389337463024351633888106D-10;
     +2.0078726603525024874472946222860416356501683654204D-11;
     +2.6961046387506385537174582238005387597979915077027D-13;
     -1.6760795833721868880525746902531471557190633459574D-14;
     -1.7118124690480244785507671262225642919352327033033D-16;
     +9.1839792232715994576195555638833867759620923522355D-18];
  fodd=zeros(size(x));
  ind=find(x<4/5);
  if ~isempty(ind)
    b = b/2;
    x2 = x(ind).^2;
    m = 17;
    tmp = b(m-1)*ones(size(x2));
    for i=m-3:-2:1
      tmp = tmp.*x2+b(i);
    end
    fodd(ind)=tmp;
  end
  ind=find(x>=4/5 & x<=10);
  if ~isempty(ind)
    nodd  = 13; 
    L = 10/sqrt(2);
    t = x(ind).^2/L^2-1;
    fodd(ind)  = chebpols(t,nodd)*codd0; 
  end

  function flog =jm1log(x)
  clog0 = [-.5001323317427930E+01;  -.7493821150778837E+01;
	   -.5102067942811160E+00;  0.2508222355964408E+01;
	   -.4427800206422171E+00;  0.3019820668279420E-01;
	   -.1067673263388026E-02;  0.2247303511677814E-04;
	   -.3065868576163888E-06;  0.2875538347259317E-08;
	   -.1937394988871976E-10;  0.9702899915931817E-13;
	   -.3712968268530488E-15;  0.1110500759672965E-17;
	   -.2625105260334293E-20];
  a=[-2.0000000000000000000000000000000000000000000000000D0;
     +1.0000000000000000000000000000000000000000000000000D0;
     -4.1666666666666666666666666666666666666666666666667D-2;
     +4.6296296296296296296296296296296296296296296296296D-4;
     -2.0667989417989417989417989417989417989417989417989D-6;
     +4.5928865373309817754262198706643151087595532039976D-9;
     -5.7990991632966941608916917558892867534842843484819D-12;
     +4.5518831737022717118459118963024228834256549046169D-15;
     -2.3707724863032665165864124459908452517841952628213D-18];
  flog=zeros(size(x));
  ind=find(x<4/5);
  if ~isempty(ind)
    a = a/2;
    x2 = x(ind).^2;
    m2 = 9;
    tmp = a(m2)*ones(size(x2));
    for i=m2-1:-1:1
      tmp=tmp.*x2+a(i);
    end
    flog(ind)=tmp;
  end
  ind=find(x>=4/5 & x<=10);
  if ~isempty(ind)
    nlog  = 14;
    L = 10/sqrt(2);
    t = x(ind).^2/L^2-1;
    flog(ind)  = chebpols(t,nlog)*clog0; 
  end

  function pols=chebpols(x,n)
  pols=zeros(length(x),n+1);
  pols(:,1)=1;
  if n == 0
    return;
  end
  pols(:,2)=x;
  for k=1:n-1
    pols(:,k+2)=2*x.*pols(:,k+1)-pols(:,k);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  functions for adaptive refinement for self and neighbor interactions
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [sinter,sinterdiff,npan]=rightdivide(t1,t2,hmax)
  dt=t2-t1;
  if  dt>hmax
    nlevel=floor(log2(dt/hmax));
    sinter=zeros(nlevel+2,1);
    sinter(1:2)=[t1;t2];
    for i=1:nlevel
      sinter(3:end)=sinter(2:end-1);
      sinter(2)=(sinter(1)+sinter(2))/2;
    end
  else
    sinter=[t1;t2];
  end
  sinterdiff=diff(sinter);
  npan=length(sinterdiff);

  function [sinter,sinterdiff,npan]=leftdivide(t1,t2,hmax)
  dt=t2-t1;
  if dt>hmax
    nlevel=floor(log2(dt/hmax));
    sinter=zeros(nlevel+2,1);
    sinter(end-1:end)=[t1;t2];
    for i=1:nlevel
      sinter(1:end-2)=sinter(2:end-1);
      sinter(end-1)=(sinter(end-2)+sinter(end))/2;
    end
  else
    sinter=[t1;t2];
  end
  sinterdiff=diff(sinter);
  npan=length(sinterdiff);  
    
  function [sinter,sinterdiff,npan,tind]=selfdivide(t1,t2,targ,hmax)
  dt=t2-t1;
  dd=dt;
  if dt>hmax
    nlevel=ceil(log2(dt/hmax));
    h=dt/2^nlevel;
    npanmax=1000;
    tt=zeros(npanmax,1);
    tt(1:2)=[t1;t2];
    tl=t1;
    tr=t2;
    ind=2;
    ilength=2;
% bisection for adaptive refinement
    while dd>h
      tc=(tl+tr)/2;
      ilength=ilength+1;
      tt(ind+1:ilength)=tt(ind:ilength-1);
      tt(ind)=tc;
      if targ>tl && targ<=tc
        tr=tc;
      else
        tl=tc;
        ind=ind+1;
      end
      dd=tr-tl;
    end
% loop to make sure the refinement is level restricted
    maxit=100;
    npan=ilength-1;
    for k=1:maxit
      npanold=npan;
      for j=1:npanold
        lj=tt(j+1)-tt(j);
        if j>1
          ll=tt(j)-tt(j-1);
        else
          ll=dt;
        end
        if j<npanold
          lr=tt(j+2)-tt(j+1);
        else
          lr=dt;
        end
        if lj>2.05*ll || lj>2.05*lr
          npan=npan+1;
          tc=(tt(j)+tt(j+1))/2;
          tt(j+2:npan+1)=tt(j+1:npan);
          tt(j+1)=tc;
        end
      end
    end
    sinter=tt(1:npan+1);
    tind=0;
    for i=1:npan
      if targ>sinter(i)
        tind=tind+1;
      end
    end
    sinterdiff=diff(sinter);
  else
    npan=1;
    tind=1;
    sinter=[t1;t2];
    sinterdiff=diff(sinter);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Legendre polynomial subroutines
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  function pols = legepols(x,n)
  % return values of legendre polynomial p_0(x), p_1(x), ..., p_n(x)
  m=length(x);
  pols=ones(m,n+1);
  if n == 0 
    return
  end
  pols(:,2)=x;
  for k=1:n-1
    pols(:,k+2)=( (2*k+1)*x.*pols(:,k+1) - k*pols(:,k) )/(k+1);
  end

  function [nodes,whts,u,v] = legeexps(n)
% *** nodes and whts are Legendre nodes of weights of length n;
% *** u is n x n matrix converting function values at n Legendre
% *** nodes to Legendre expansion coefficients.
% *** v is the inverse of u
  [nodes, whts]=legewhts(n);
  if nargout < 3
    return
  end
  u=zeros(n);
  for i=1:n
    u(:,i)=legepols(nodes(i),n-1);
  end
  v=u';
  for i=1:n
    d=(2*i-1)/2;
    for j=1:n
      u(i,j)=v(j,i)*whts(j)*d;
    end
  end

  function [ts, whts]=legewhts(n)
  eps=1d-15;
  h=pi/(2*n);
  ts=cos(h*(1:2:2*n-1));
  ts(floor(n/2)+1)=0;
  for i=1:floor(n/2)
    xk=ts(i);
    ifout=0;
    for k=1:10
      [pol,der,~]=legepol_sum(xk,n);
      delta=-pol/der;
      xk=xk+delta;
      if abs(delta) < eps
        ifout=ifout+1;
      end
      if ifout == 3
        break;
      end
    end
    ts(i)=xk;
    ts(n-i+1)=-xk;
  end
  ts=sort(ts);
  ts=ts';
  if nargout == 1
    return
  end
  whts=zeros(n,1);
  for i=1:floor((n+1)/2)
    [~,~,sum]=legepol_sum(ts(i),n);
    whts(i)=1/sum;
    whts(n-i+1)=whts(i);
  end

  function [pol,der,sum]=legepol_sum(x,n)
  sum=0;
  pol=1;
  der=0;
  sum=sum+pol^2/2;
  if n==0
    return;
  end
  pol=x;
  der=1;
  sum=sum+pol^2*1.5;
  if n==1
    return
  end
  sum=0;
  pkm1=1;
  pk=x;
  sum=sum+pkm1^2/2;
  sum=sum+pk^2*1.5;
  pk=1;
  pkp1=x;
  for k=1:n-1
    pkm1=pk;
    pk=pkp1;
    pkp1=( (2*k+1)*x*pk-k*pkm1)/(k+1);
    sum=sum+pkp1^2*(k+1.5);
  end
  pol=pkp1;
  der=n*(x*pkp1-pk)/(x^2-1);
