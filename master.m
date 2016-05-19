function master


Re = 64;
cfl = 50;
imax = 100000;
tol = 1e-12;
k=4;

x0 = -4;
xL =  4;
n = 17;
x = linspace(x0,xL,n)';
delx = (xL-x0)/(n-1);
for i = 2:n-1
    x(i) = x(i) + 0.2*delx*(2*rand(1)-1);
end
dx = x(2:end) - x(1:end-1);
xc = (x(2:end)+x(1:end-1))/2;

fprintf('Normalized minimum area: %3.1f\n',min(dx)/delx*100);
fprintf('Normalized maximum area: %3.1f\n',max(dx)/delx*100);
figure(1)
plot(xc,dx,'k-o')
title('grid spacing')

if any(dx<eps)
    fprintf('Negative Area!!!!! Ahhhhhhhh!!!!\n')
    pause
end


lref = 8;%xL-x0;
uref = 2;
nu = uref*lref/Re; %Re = u*l/nu

%Cosine exact
a=4; b=3.5; c=1; d=2;
fex = @(x) a + b*cos(c*x-d*pi);
dfex = @(x) -c*b*sin(c*x-d*pi);

% Burger's exact
% c1=-uref; c2=Re/(2*lref);
% fex = @(x) c1*tanh( c2 * x );
% dfex = @(x) -c2/c1*(fex(x).^2-c1.^2);

source = sourceterm(x,nu,fex,dfex);

uexact = fexact(x, fex);

usoln = burgers(x, nu, source, cfl, tol, imax, uexact, fex);

de = usoln(2:end-1) - uexact;
l2de = sqrt( de'*de/numel(de) );
fprintf('%23.15e\n',l2de)

te_exact = exact_truncation_error(x, nu, source, uexact, fex);
de_ete_test = burgers_ete(x, usoln, te_exact, nu, 1e9, tol, 1000);
fprintf('ETE passed with exact te? %5.3f %%\n',max( abs(de_ete_test(2:end-1)-de)./abs(de) ) );

te = te_estimate(usoln, source, x, k, nu);
desoln = burgers_ete(x, usoln, te, nu, 1e9, tol, 1000);

figure(2)
subplot(3,1,1)
plot(xc,uexact,'k-o',xc,usoln(2:end-1),'r-o')
title('Exact')

subplot(3,1,2)
plot(xc,te_exact,'k-o',xc,te,'r-*')
title('TE estimate')

subplot(3,1,3)
plot(xc,de,'k-o',xc,de_ete_test(2:end-1),'k-*',xc,desoln(2:end-1),'r-o')
title('DE estimate')

tees_err(:,1)=te;
de_diff = 1;
n = 2;

% while n<=10;
% 
% 
%     % Estimate dte/du
%     utest = usoln(2:end-1)-desoln(2:end-1,n-1);
% %     te3 = exact_truncation_error(x, nu, source, utest, fex);
%     for i = 1:length(uexact)
%         u1 = utest; u1(i)=u1(i)+1e-7;
%         u2 = utest; u2(i)=u2(i)-1e-7;
%         %   u1 = uexact; u1(i)=u1(i)+1e-7;
%         %   u2 = uexact; u2(i)=u2(i)-1e-7;
%         %   te1 = te_estimate(u1, x, k, nu);
%         %   te2 = te_estimate(u2, x, k, nu);
%         te1 = exact_truncation_error(x, nu, source, u1, fex);
%         te2 = exact_truncation_error(x, nu, source, u2, fex);
%         dtedu(:,i) = -(  te1 - te2 )'/2e-7;
%     end
% 
%     % Estimate dte/du
% 
% %     utest = usoln;%-desoln;
% %     for i = 2:length(utest)-1
% %       u1 = utest; u1(i)=u1(i)+1e-7;
% %       u2 = utest; u2(i)=u2(i)-1e-7;
% %       te1 = te_estimate(u1, source, x, k, nu);
% %       te2 = te_estimate(u2, source, x, k, nu);
% % 
% %       dtedu(:,i-1) = ( te1 - te2 )'/2e-7;    
% %     end
% %     utest = utest(2:end-1);
%     
% 
%     tees_err(:,n) = dtedu*desoln(2:end-1,n-1);
% %     tees_err(:,n) = te3;
%     [desoln1, conv] = burgers_ete(x, usoln, tees_err(:,n), nu, 1e9, tol, 1000);
%     desoln(:,n)=desoln1;
%     
%     de_diff = max(  abs(desoln(:,n)-desoln(:,n-1)) );
%     fprintf('DE tolerance %8.4e\n',de_diff)
% 
%     figure(3)
%     subplot(3,1,1)
%     plot(xc,uexact,'k-o',xc,usoln(2:end-1),'r-o',xc,utest,'g-v')
%     title('Solution')
% 
%     subplot(3,1,2)
%     plot(xc,te_exact,'k-o',xc,te,'r-o',xc,tees_err(:,n),'g-v')
%     title('TE estimate')
% 
%     subplot(3,1,3)
%     plot(xc,de,'k-o',xc,desoln(2:end-1,1),'r-o',xc,desoln(2:end-1,n),'g-v')
%     title('DE estimate')
% 
%     n = n + 1;
% end

% tees_err = dtedu*desoln(2:end-1);
% tees_err = dtedu*de;

% figure(4)
% subplot(3,1,1)
% plot(xc,uexact,'k-o',xc,utest,'g-v')
% title('Solution')
% 
% subplot(3,1,2)
% plot(xc,te_exact,'k-o',xc,tees_err(:,n-1),'g-v')
% title('TE estimate')
% 
% subplot(3,1,3)
% plot(xc,de,'k-o',xc,desoln(2:end-1,n-1),'g-v')
% title('DE estimate')

%% This is a test sample u = uh + de
% Random sampleing of u with a variance of de to compute te.
% figure(2)
% hold off
% plot(xc,te_exact,'k-o',xc,te,'r-o')
% hold on
% title('TE estimate')
% 
% 
% teuex = te_estimate(usoln - desoln, x, k, nu);
% 
% 
% 
% 

% Estimate dte/du
utest = usoln(2:end-1)-desoln(2:end-1,n-1);
%     te3 = exact_truncation_error(x, nu, source, utest, fex);
for i = 1:length(uexact)
    u1 = utest; u1(i)=u1(i)+1e-7;
    u2 = utest; u2(i)=u2(i)-1e-7;
    %   u1 = uexact; u1(i)=u1(i)+1e-7;
    %   u2 = uexact; u2(i)=u2(i)-1e-7;
    %   te1 = te_estimate(u1, x, k, nu);
    %   te2 = te_estimate(u2, x, k, nu);
    te1 = exact_truncation_error(x, nu, source, u1, fex);
    te2 = exact_truncation_error(x, nu, source, u2, fex);
    dtedu(:,i) = -(  te1 - te2 )'/2e-7;
end

nsamples = 300;
cnt = 2;
detest(:,1) = desoln;
de_mean = mean(detest,2);
u_ex_est = usoln - desoln;

% Randomly sample de variations to determine how sensitive the solution is
% to de
% for i = 1:nsamples
% %         de_eff =  desoln.*(2*rand(size(desoln))-1);
%         de_eff =  desoln.*(randn(size(desoln)));
% %         de_eff =  desoln.*(randn(size(desoln)));
%         utest = usoln - de_eff;
% %         tetest(:,i) = te_estimate(utest, source, x, k, nu);
%         tetest(:,i) = exact_truncation_error(x, nu, source, utest(2:end-1), fex);
%     
%     [detest1, conv] = burgers_ete(x, usoln, tetest(:,i), nu, 1e9, 1e-9, 100);
%     if (conv==1); detest2(:,cnt) = detest1; cnt = cnt + 1; end;
% 
%     figure(4)
%     plot(xc, de, 'k-o', xc, desoln(2:end-1), 'g-*',xc, detest1(2:end-1), 'b-*',xc, de_eff(2:end-1), 'b-v',xc, detest2(2:end-1,cnt-1), 'b-o', xc, de_mean(2:end-1),'r-o')
% %     hold on
%     de_mean = mean(detest2,2);
%     de_std = std(detest2,0,2);
%     
%     te_mean = mean(tetest,2);
%     te_std = std(tetest,0,2);
% end


nsamples = 100;
cnt = 2;
detest(:,1) = desoln;
de_mean = mean(detest,2);
u_ex_est = usoln - desoln;

    tetest1 = te_estimate(usoln, source, x, k, nu);
    tetest2 = te_estimate(usoln, source, x, k+1, nu);
    dte = abs(tetest1 - tetest2);
    
%     tetest3 = exact_truncation_error(x, nu, source, usoln(2:end-1)-desoln(2:end-1), fex);
for i = 1:nsamples
%         de_mean = mean(detest,2);
%         de_std = std(detest(:,1:50),0,2);
%         de_pur = de_std.*(2*rand(size(desoln))-1);
%         de_pur = desoln.*randn(size(desoln));
        
%         de_eff = desoln - de_pur;
%         uex_est = usoln(2:end-1) - de_eff(2:end-1);
%         tetest(:,i) = exact_truncation_error(x, nu, source, uex_est, fex);
%         ueff = usoln + de_pur;
%         tetest(:,i) = te_estimate(ueff, source, x, k, nu);



    tetest(:,i) = te + dte.*randn(size(te));
%     tetest(:,i) = te + dte.*(2*rand(size(te)) - 1);
%     tetest(:,i) = te + dte.*randn(1);
    
    [detest1, conv] = burgers_ete(x, usoln, tetest(:,i), nu, 1e9, 1e-9, 100);
    if (conv==1); detest(:,cnt) = detest1; cnt = cnt + 1; end;
    
    
%     de_mean = mean(detest,2);
    figure(4)
    plot(xc, de, 'k-o', xc, desoln(2:end-1), 'g-*', xc, de_mean(2:end-1),'r-o')
%     hold on
end
  te_mean = mean(tetest,2);
  te_std = std(tetest,0,2);
  de_mean = mean(detest,2);
  de_std = std(detest,0,2);
  
  [de_temean, conv] = burgers_ete(x, usoln, te_mean, nu, 1e9, 1e-9, 100);
  [de_temean1, conv] = burgers_ete(x, usoln, te_mean-te_std, nu, 1e9, 1e-9, 100);
  [de_temean2, conv] = burgers_ete(x, usoln, te_mean+te_std, nu, 1e9, 1e-9, 100);
  
  
  [de_te4, conv] = burgers_ete(x, usoln, tetest1, nu, 1e9, 1e-9, 100);
  [de_te6, conv] = burgers_ete(x, usoln, tetest2, nu, 1e9, 1e-9, 100);
  
  figure(2)
  subplot(2,1,1)
%   for i = 51:nsamples
%       if (i==1); hold off; end;
%       plot(xc, tetest(:,i),'g-*')
%       hold on
%   end
  plot(xc,te_exact,'k-o',xc,te,'r-o',xc, te_mean,'r-*')
  hold on
  errorbar(xc, te_mean, 2*te_std)
%   errorbar(xc, te, 2*te_std)
  title('TE estimate')

  subplot(2,1,2)
%   for i = 51:size(detest,2)
%       if (i==1); hold off; end;
%       plot(xc, detest(2:end-1,i),'g-*')
%       hold on
%   end
  plot(xc,de,'k-o',xc,de_te4(2:end-1),'g-o',xc,de_te6(2:end-1),'g--*',xc, de_mean(2:end-1),'r-*')
  hold on
  errorbar(xc, de_mean(2:end-1), 2*de_std(2:end-1))
%   errorbar(xc, de_te4(2:end-1), de_te4(2:end-1)-de_te6(2:end-1))
%   errorbar(xc, de_temean(2:end-1), de_temean(2:end-1)-de_temean1(2:end-1), de_temean(2:end-1)+de_temean1(2:end-1) )
%   axis([min(xc),max(xc),4*min(de),4*max(de)])
  title('DE estimate')

end

function [usoln] = burgers(x, nu, source, cfl, tol, imax, uinit, fex)

x0 = min(x);
xL = max(x);
n = numel(x);
x = reshape(x,[n,1]);


%% Input parameters
% x0 = -4;
% xL =  4; 
% n = 65;
% x = linspace(x0,xL,n);
% Re = 16;
% cfl = 0.1;
% imax = 200000;
% tol = 1e-13;



% Add ghost cells
x1 = x0 - (x(2)-x(1));
xnp1 = xL + (x(end)-x(end-1));
x = [x1; x; xnp1];

% Add ghost cells
u1 = fexact(x(1:2), fex);
un = fexact(x(end-1:end), fex);
u = [u1; uinit; un];

xc = (x(2:end)+x(1:end-1))/2;
dx = x(2:end)-x(1:end-1);



%Set up dt array
dt_arry = diag(u);
I = find(dt_arry~=0); dt_arry(I) = 0; I=I(2:end-1);
dt_arry(I) = 1;
dt_arry = sparse(dt_arry);
Itime = I;

i=1;
l2resid=1;
while i<imax && l2resid>tol
    r = resid(x,u,nu)+source;
    l2resid = sqrt( r'*r/(numel(r)) );
    l2resid_vec(i) = l2resid;
    if (mod(i,1000)==0); fprintf('%5.0f  %12.4e\n', i, l2resid); end;

    %apply bc
    b = [-(u(1)+u(2)+4*tanh(-8)); r; -(u(end)+u(end-1)+4*tanh(8))];
    b = [0; r; 0];

    
    % Compute time step
    dt = min( cfl*dx(2:end-1).^2 ./ (2*nu + abs( dx(2:end-1).*u(2:end-1) ) ) );
    
% Implicit
    df = flux_jacobian(x,u,nu);
    
    df(Itime) = df(Itime) + dx(2:end-1)/dt;
    du = df\b;
% Explicit
%     du = dt*b;
    
    unp1 = u + du;
    
    u = unp1;
    i = i + 1;
end
fprintf('%5.0f  %12.4e\n', i, l2resid);

usoln = u;


end

function [desoln, conv] = burgers_ete(x, u, te, nu, cfl, tol, imax)
% (u(n+1)-u(n))/dt - df/du*(u(n+1)-u(n)) = 0


x0 = min(x);
xL = max(x);
n = numel(x);
x = reshape(x,[n,1]);

% Add ghost cells
x1 = x0 - (x(2)-x(1));
xnp1 = xL + (x(end)-x(end-1));
x = [x1; x; xnp1];

% Initialize DE
de = zeros(size(x,1)-1,1);

xc = (x(2:end)+x(1:end-1))/2;
dx = x(2:end)-x(1:end-1);

% Jacobian
dR_du = flux_jacobian(x,u,nu);

% Set time step. It's constant...
dt = min( cfl*dx(2:end-1).^2 ./ (2*nu + abs( dx(2:end-1).*u(2:end-1) ) ) );

%Set up dt array
dt_arry = diag(u);
I = find(dt_arry~=0); dt_arry(I) = 0; I=I(2:end-1);
dt_arry(I) = 1;
dt_arry = sparse(dt_arry);
Itime = I;
dt_arry(Itime) = 1/dt;

i=1;
l2resid=1;
while i<imax && l2resid>tol
    r = te - dR_du(2:end-1,2:end-1)*de(2:end-1);

    l2resid = sqrt( r'*r/(numel(r)) );
    l2resid_vec(i) = l2resid;
    if (mod(i,100)==0); fprintf('%5.0f  %12.4e\n', i, l2resid); end;

    %apply bc
    b = [0; r; 0];

    df = dR_du + dt_arry;
    dde = df\b;

    denp1 = de + dde;

    de = denp1;
    i = i + 1;
end
% fprintf('%5.0f  %12.4e\n', i, l2resid);
if l2resid>tol; conv = 0; else conv = 1; end

% desoln = de(2:end-1);
desoln = de;

end

function f = flux(x,u,nu)
  dx = x(2:end)-x(1:end-1);
  i = 1:length(u)-1;
  f(i,1) = (u(i+1).^2 + u(i).^2)/4 - nu* (u(i+1)-u(i))./( 1/2*(dx(i+1)+dx(i)) );
end

function f = analytic_flux(x,u,du,nu)
  f = u.^2/2 - nu*du;
end

function r = resid(x,u,nu)
  dx = x(2:end)-x(1:end-1);
  f = flux(x, u, nu);
  r = -( f( 2:end ) - f( 1:end-1 ) );
end

function df = flux_jacobian(x,u,nu)
% compute the flux jacobian for
% f(i+1/2) = [u(i+1).^2 + u(i)^2 ]/4 - nu*2*(u(i+1)-u(i))/( dx(i+1)+dx(i) );
% f(i-1/2) = [u(i).^2 + u(i-1)^2 ]/4 - nu*2*(u(i)-u(i-1))/( dx(i)+dx(i-1) );
% rhs = (f(i+1/2) - f(i-1/2))/dx(i)
% rhs = (u(i+1).^2 - u(i-1).^2)/4 - nu*2*(u(i+1)-u(i))/( dx(i+1)+dx(i) ) + nu*2*(u(i)-u(i-1))/( dx(i)+dx(i-1) )
  i = 2:length(u)-1;
  dx = x(2:end)-x(1:end-1);
  dfdu_ip1(i) = u(i+1)/2 - nu*2./(dx(i+1)+dx(i)) ; dfdu_ip1(1) = 0; 
  dfdu_im1(i-1) = u(i-1)/2 + nu*2./( dx(i)+dx(i-1) ); dfdu_im1(length(u)-1)=0;
  dfdu_i(i)   = nu*2./(dx(i+1)+dx(i)) + nu*2./( dx(i) + dx(i-1) ); dfdu_i(1) = 1; dfdu_i(length(u))=1;
   
  df = (diag(dfdu_i) + diag(dfdu_ip1,1) - diag(dfdu_im1,-1));
end

function fex = fexact(x, fex_fun)
  dx = x(2:end) - x(1:end-1);
  for i = 1:length(dx)
    fex(i,1) = integral(fex_fun, x(i), x(i+1))/(x(i+1)-x(i));
  end
end

function s = sourceterm(x, nu, fex_fun, dfex_fun)
  dx = x(2:end) - x(1:end-1);
  u = fex_fun(x);
  du = dfex_fun(x);

  f = analytic_flux( x, u, du, nu );
  s = (f(2:end)-f(1:end-1));
end

function te_exact = exact_truncation_error(x, nu, source, u, fex)
    x0 = min(x);
    xL = max(x);
    n = numel(x);
    x = reshape(x,[n,1]);

    % Add ghost cells
    x1 = x0 - (x(2)-x(1));
    xnp1 = xL + (x(end)-x(end-1));
    x = [x1; x; xnp1];

    % Add ghost cells
    u1 = fexact(x(1:2), fex);
    un = fexact(x(end-1:end), fex);
    u = [u1; u; un];

    te_exact = resid(x,u,nu)+source;
end

function te = te_estimate(u, source, x, k, nu)
    x0 = min(x);
    xL = max(x);
    n = numel(x);
    x = reshape(x,[n,1]);

    % Add ghost cells
    x1 = x0 - (x(2)-x(1));
    xnp1 = xL + (x(end)-x(end-1));
    x = [x1; x; xnp1];

  p = 0:k;
  dx = (x(2:end)-x(1:end-1));

  fun = @(x,a) sum(a.*x.^p);
  dfun = @(x,a) sum(p(2:end).*a(2:end).*x.^(p(2:end)-1));
  flux = @(x,a) fun(x,a).^2/2 - nu*dfun(x,a);

  for i = 2:length(u)-1
     il = i-floor( (k+1)/2 );
     ih = il + k;
     
     if (il<1); il = 1; ih = il+k; end
     if (ih>numel(u)); ih = numel(u); il = ih - k; end
          
     b = u(il:ih);
     
     for j = 1:k+1
        A(j,:) = (x(il+j).^(p+1)./(p+1) - x(il+j-1).^(p+1)./(p+1))./(dx(il+j-1));
     end
     
     a = A\b;
     te(i-1,1) = flux(x(i+1),a') - flux(x(i),a');
  end

  te = te - source;

end
