function te_burgers_analysis
%% Burgers TE analysis
% F = u^2/2 - nu*du/dx

re=32; l=8; u=2; nu = l*u/re;
c1=-u; c2 = re/(2*l);

u = @(x) c1*tanh(c2*x);
ux = @(x) -c2/c1*(u(x).^2 - c1.^2);
uxx = @(x) -c2/c1*(2*u(x).*ux(x));
uxxx = @(x) -2*c2/c1*(ux(x).^2 + u(x).*uxx(x));
uxxxx = @(x) -2*c2/c1*(3*ux(x).*uxx(x)+u(x)*uxxx(x));
uxxxxx = @(x) -2*c2/c1*(3*(uxx(x).^2+ux(x).*uxxx(x))+ux(x).*uxxx(x)+u(x).*uxxxx(x));
uxxxxxx = @(x) -2*c2/c1*( 6*uxx(x)*uxxx(x) + 5*ux(x).*uxxxx(x) + 4*uxx(x).*uxxx(x) + u(x).*uxxxxx(x) );
ubar = @(x1,x2) integral(u, x1, x2, 'AbsTol',1e-10)/(x2-x1);

%Test analytic derivatives (p==2)
x0=1; dx=0.5; n = 10;
ux_p = DerivativeTest(u,ux,x0,dx,n);
uxx_p = DerivativeTest(ux,uxx,x0,dx,n);
uxxx_p = DerivativeTest(uxx,uxxx,x0,dx,n);
uxxxx_p = DerivativeTest(uxxx,uxxxx,x0,dx,n);
uxxxxx_p = DerivativeTest(uxxxx,uxxxxx,x0,dx,n);
uxxxxxx_p = DerivativeTest(uxxxxx,uxxxxxx,x0,dx,n);

% Taylor Series expansion test
% F = @(x) u(x);
% Fh = @(x,dx) u(x+dx/2);
% TE = @(x,dx) dx/2*ux(x) + dx^2/8*uxx(x) + dx^3/48*uxxx(x) + dx^4/384*uxxxx(x);
% pout=OrderTest(x0, dx, F, Fh, TE, 10)



% Analytic TE for u^2/2
%      6               6               6               6     2     4             4             4    2     2           2   2    2
%    dx  u uxxxxxx   dx  ux uxxxxx   dx  uxx uxxxx   dx  uxxx    dx  u uxxxx   dx  ux uxxx   dx  uxx    dx  u uxx   dx  ux    u
%    ------------- + ------------- + ------------- + --------- + ----------- + ----------- + -------- + --------- + ------- + --
%         46080            7680            3072          4608         384            96         128         8          8       2


% Centered flux average
fprintf('Discrete Face Flux\n')
F = @(x) u(x).^2/2;
Fh = @(x,dx) (u(x+dx/2).^2+u(x-dx/2).^2)/4;
TE2 = @(x,dx) dx.^2*( u(x).*uxx(x)/8 + ux(x).^2/8); % + O(dx^4)
TE4 = @(x,dx) dx.^4*( uxx(x).^2/128 + ux(x).*uxxx(x)/96 + u(x).*uxxxx(x)/384 ); % + O(dx^6)
TE = @(x,dx) TE2(x,dx) + TE4(x,dx);

pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('Central average O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) , 7); fprintf('Central average O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) +TE4(x,dx) , 5); fprintf('Central average O(dx^6): %4.2f\n\n',pout);

fprintf('Discrete Face Flux using the FV solution approximation |  o--x--o  |\n')
F = @(x) u(x).^2/2 - nu*ux(x);
Fh = @(x,dx) ( ubar(x,x+dx).^2+ubar(x-dx,x).^2)/4 - nu*(ubar(x,x+dx) - ubar(x-dx,x))/(dx);
TE2 = @(x,dx) dx.^2*( u(x).*uxx(x)/6 + ux(x).^2/8 - nu*uxxx(x)/12); % + O(dx^4)
TE4 = @(x,dx) dx.^4*( uxx(x).^2/72 + ux(x).*uxxx(x)/48 + u(x).*uxxxx(x)/120 -nu*uxxxxx(x)/360); % + O(dx^6)
TE = @(x,dx) TE2(x,dx) + TE4(x,dx);

pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('Central average O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) , 7); fprintf('Central average O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) +TE4(x,dx) , 5); fprintf('Central average O(dx^6): %4.2f\n\n',pout);


% Forward flux average
fprintf('Discrete Face Flux approximation |  x--|--o  |\n')
F = @(x) u(x).^2/2;
Fh = @(x,dx) (u(x+dx).^2+u(x).^2)/4;
TE1 = @(x,dx) dx*( u(x).*ux(x)/2 ); % + O(dx^2)
TE2 = @(x,dx) dx.^2*( ux(x).^2/4 + u(x).*uxx(x)/4 ); % + O(dx^3)
TE3 = @(x,dx) dx.^3*( ux(x).*uxx(x)/4 + u(x).*uxxx(x)/12); % + O(dx^4)
TE4 = @(x,dx) dx.^4*( u(x).*uxxxx(x)/48 + ux(x).*uxxx(x)/12 + uxx(x).^2/16); % + O(dx^5)

TE = @(x,dx) TE1(x,dx) +TE2(x,dx) +TE3(x,dx) +TE4(x,dx);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('Forward average O(dx^1): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) , 10); fprintf('Forward average O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) +TE2(x,dx) , 10); fprintf('Forward average O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) +TE2(x,dx) +TE3(x,dx)  , 9); fprintf('Forward average O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) +TE2(x,dx) +TE3(x,dx) +TE4(x,dx) , 7); fprintf('Forward average O(dx^5): %4.2f\n\n',pout);

% Forward flux average
fprintf('Discrete Face Flux using the FV solution approximation |  x--|--o  |\n')
F = @(x) u(x).^2/2 - nu*ux(x);
Fh = @(x,dx) ( ubar(x+dx/2,x+3/2*dx).^2+ubar(x-dx/2,x+dx/2).^2)/4 - nu*(ubar(x+dx/2,x+3/2*dx)-ubar(x-dx/2,x+dx/2))/dx;
TE1 = @(x,dx) dx*( u(x).*ux(x)/2 -nu*uxx(x)/2); % + O(dx^2)
TE2 = @(x,dx) dx.^2*( ux(x).^2/4 + u(x).*uxx(x)*7/24 -nu*uxxx(x)*5/24); % + O(dx^3)
TE3 = @(x,dx) dx.^3*( ux(x).*uxx(x)*13/48 + u(x).*uxxx(x)*5/48 - nu*uxxxx(x)*1/16); % + O(dx^4)
TE4 = @(x,dx) dx.^4*( u(x).*uxxxx(x)*61/1920 + ux(x).*uxxx(x)*5/48 + uxx(x).^2*85/1152 - nu*uxxxxx(x)*91/5760); % + O(dx^5)

TE = @(x,dx) TE1(x,dx) +TE2(x,dx) +TE3(x,dx) +TE4(x,dx);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('Forward average O(dx^1): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) , 10); fprintf('Forward average O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) +TE2(x,dx) , 10); fprintf('Forward average O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) +TE2(x,dx) +TE3(x,dx)  , 9); fprintf('Forward average O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) +TE2(x,dx) +TE3(x,dx) +TE4(x,dx) , 7); fprintf('Forward average O(dx^5): %4.2f\n\n',pout);


% Backward flux average
fprintf('Discrete Face Flux approximation |  o--|--x  |     |\n')
F = @(x) u(x).^2/2;
Fh = @(x,dx) (u(x-dx).^2+u(x).^2)/4;
TE1 = @(x,dx) -dx*( u(x).*ux(x)/2 ); % + O(dx^2)
TE2 = @(x,dx) dx.^2*( ux(x).^2/4 + u(x).*uxx(x)/4 ); % + O(dx^3)
TE3 = @(x,dx) -dx.^3*( ux(x).*uxx(x)/4 + u(x).*uxxx(x)/12); % + O(dx^4)
TE4 = @(x,dx) dx.^4*( u(x).*uxxxx(x)/48 + ux(x).*uxxx(x)/12 + uxx(x).^2/16); % + O(dx^5)

TE = @(x,dx) TE1(x,dx) +TE2(x,dx) +TE3(x,dx) +TE4(x,dx);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('Backward average O(dx^1): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) , 10); fprintf('Backward average O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) +TE2(x,dx) , 10); fprintf('Backward average O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) +TE2(x,dx) +TE3(x,dx)  , 9); fprintf('Backward average O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) +TE2(x,dx) +TE3(x,dx) +TE4(x,dx) , 7); fprintf('Backward average O(dx^5): %4.2f\n\n',pout);

fprintf('Discrete Face Flux using the FV solution approximation |  o--|--x  |     |\n')
F = @(x) u(x).^2/2-nu*ux(x);
Fh = @(x,dx) ( ubar(x-dx*3/2,x-1/2*dx).^2+ubar(x-dx/2,x+dx/2).^2)/4 - nu*(ubar(x-dx/2,x+dx/2)-ubar(x-dx*3/2,x-1/2*dx))/dx;
TE1 = @(x,dx) -dx*( u(x).*ux(x)/2 -nu*uxx(x)/2); % + O(dx^2) <---- TE1=0 cause it's the continuous operator (by definition)
TE2 = @(x,dx) dx.^2*( ux(x).^2/4 + u(x).*uxx(x)*7/24 -nu*uxxx(x)*5/24); % + O(dx^3)
TE3 = @(x,dx) -dx.^3*( ux(x).*uxx(x)*13/48 + u(x).*uxxx(x)*5/48 - nu*uxxxx(x)*1/16); % + O(dx^4)
TE4 = @(x,dx) dx.^4*( u(x).*uxxxx(x)*61/1920 + ux(x).*uxxx(x)*5/48 + uxx(x).^2*85/1152 - nu*uxxxxx(x)*91/5760); % + O(dx^5)

TE = @(x,dx) TE1(x,dx) +TE2(x,dx) +TE3(x,dx) +TE4(x,dx);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('Backward average O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) , 10); fprintf('Backward average O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) +TE2(x,dx) , 10); fprintf('Backward average O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) +TE2(x,dx) +TE3(x,dx)  , 9); fprintf('Backward average O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) +TE2(x,dx) +TE3(x,dx) +TE4(x,dx) , 7); fprintf('Backward average O(dx^5): %4.2f\n\n',pout);


% Centered flux difference
fprintf('Discrete cell residual approximation |--o--|  x  |--o--|\n')
F = @(x) u(x).*ux(x);
Fh = @(x,dx) (u(x+dx).^2-u(x-dx).^2)/(4*dx);
TE2 = @(x,dx) dx.^2*( u(x).*uxxx(x)/6 + ux(x).*uxx(x)/2); % + O(dx^4)
TE4 = @(x,dx) dx.^4*( ux(x).*uxxxx(x)/24 + uxx(x).*uxxx(x)/12 +u(x).*uxxxxx(x)/120 ); % + O(dx^6)
TE = @(x,dx) TE2(x,dx) + TE4(x,dx);

pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('Central flux difference O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) , 7); fprintf('Central flux difference O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) +TE4(x,dx) , 5); fprintf('Central flux difference O(dx^6): %4.2f\n\n',pout);



fprintf('Discrete cell residual the FV solution approximation |--o--|  x  |--o--|\n')
F = @(x) u(x).*ux(x) - nu*uxx(x);
Fh = @(x,dx) ( ubar(x+dx/2,x+dx*3/2).^2 - ubar(x-dx*3/2,x-1/2*dx).^2)/(4*dx) - nu*(ubar(x+dx/2,x+dx*3/2)-2*ubar(x-dx*1/2,x+1/2*dx) + ubar(x-dx*3/2,x-1/2*dx))/dx^2;
TE2 = @(x,dx) dx.^2*( u(x).*uxxx(x)*5/24 + ux(x).*uxx(x)*13/24 - nu*uxxxx(x)/8); % + O(dx^4)
TE4 = @(x,dx) dx.^4*( ux(x).*uxxxx(x)*121/1920 + uxx(x).*uxxx(x)*65/576 + u(x)*uxxxxx(x)*91/5760 - nu*uxxxxxx(x)*13/1920 ); % + O(dx^6)
TE = @(x,dx) TE2(x,dx) + TE4(x,dx);

pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('Central flux difference O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) , 7); fprintf('Central flux difference O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) +TE4(x,dx) , 5); fprintf('Central flux difference O(dx^6): %4.2f\n\n',pout);

% fprintf('Discrete cell residual the FV solution approximation (Dirichlet inconsistent BC: |--o--|  x  |--o--|\n')


% fprintf('u^2/2 term in the BC flux\n')
% u2A = @(x) u(x)^2/2;
% u2 = @(x,dx) u(x-dx/2)^2/2;
% TE1 = @(x,dx) -dx*(u(x)*ux(x)/2);
% TE2 = @(x,dx) dx^2*( ux(x)^2/8 + u(x)*uxx(x)/8 );
% TE3 = @(x,dx) -dx^3*( ux(x)*uxx(x)/16 +u(x)*uxxx(x)/48 );
% TE4 = @(x,dx) dx^4*( uxx(x)^2/128 + ux(x)*uxxx(x)/96 + u(x)*uxxxx(x)/384 );
% TE5 = @(x,dx) -dx^5*( uxx(x)*uxxx(x)/384 + ux(x)*uxxxx(x)/768 + u(x)*uxxxxx(x)/3840 );
% pout=OrderTest(x0, dx, u2A, u2, @(x,dx) 0 , 10); fprintf('Non linear term O(dx^1): %4.2f\n',pout);
% pout=OrderTest(x0, dx, u2A, u2, @(x,dx) TE1(x,dx) , 9); fprintf('Non linear term O(dx^2): %4.2f\n',pout);
% pout=OrderTest(x0, dx, u2A, u2, @(x,dx) TE1(x,dx) + TE2(x,dx) , 8); fprintf('Non linear term O(dx^3): %4.2f\n',pout);
% pout=OrderTest(x0, dx, u2A, u2, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx), 7); fprintf('Non linear term O(dx^4): %4.2f\n',pout);
% pout=OrderTest(x0, dx, u2A, u2, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx), 6); fprintf('Non linear term O(dx^5): %4.2f\n',pout);
% pout=OrderTest(x0, dx, u2A, u2, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx), 5); fprintf('Non linear term O(dx^6): %4.2f\n\n',pout);


% BTE_Dirichlet_inconsistent_3rd
% BTE_Dirichlet_inconsistent_4th
% 
% BTE_Neumann_consistent_3rd
% BTE_Neumann_consistent_4th
% 
% BTE_Neumann_inconsistent_3rd
BTE_Neumann_inconsistent_4th






% FV solution average approximation (centered)

fprintf('FV approximation cell centered refernce |     |--x--|     |\n')
F = @(x) u(x);
Fh = @(x,dx) ubar(x-dx/2,x+dx/2);
TE2 = @(x,dx) uxx(x)/24*dx^2;
TE4 = @(x,dx) uxxxx(x)/1920*dx^4;
TE6 = @(x,dx) uxxxxxx(x)/322560*dx^6;

pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('FV solution average O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, TE2 , 7); fprintf('FV solution average O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) + TE4(x,dx) , 4); fprintf('FV solution average O(dx^6): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) + TE4(x,dx) + TE6(x,dx) , 4); fprintf('FV solution average O(dx^8): %4.2f\n\n',pout);

% FV solution average approximation (forward cell average at left face)
fprintf('FV approximation with a forward integral |     |     x-----|\n')
F = @(x) u(x);
Fh = @(x,dx) ubar(x,x+dx);
TE1 = @(x,dx) ux(x)/2*dx;
TE2 = @(x,dx) uxx(x)/6*dx^2;
TE3 = @(x,dx) uxxx(x)/24*dx^3;
TE4 = @(x,dx) uxxxx(x)/120*dx^4;
TE5 = @(x,dx) uxxxxx(x)/720*dx^5;
TE6 = @(x,dx) uxxxxxx(x)/5040*dx^6;

pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('FV solution average O(dx^1): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, TE1 , 10); fprintf('FV solution average O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) , 10); fprintf('FV solution average O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) , 9); fprintf('FV solution average O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) , 7); fprintf('FV solution average O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx) , 6); fprintf('FV solution average O(dx^6): %4.2f\n\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx) +TE6(x,dx) , 6); fprintf('FV solution average O(dx^7): %4.2f\n\n',pout);


% FV solution average approximation (backward cell average at right face)
fprintf('FV approximation with a backward integral |-----x     |     |\n')
F = @(x) u(x);
Fh = @(x,dx) ubar(x-dx,x);
TE1 = @(x,dx) -ux(x)/2*dx;
TE2 = @(x,dx) uxx(x)/6*dx^2;
TE3 = @(x,dx) -uxxx(x)/24*dx^3;
TE4 = @(x,dx) uxxxx(x)/120*dx^4;
TE5 = @(x,dx) -uxxxxx(x)/720*dx^5;

pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('FV solution average O(dx^1): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, TE1 , 10); fprintf('FV solution average O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) , 10); fprintf('FV solution average O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) , 9); fprintf('FV solution average O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) , 7); fprintf('FV solution average O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx) , 6); fprintf('FV solution average O(dx^6): %4.2f\n\n',pout);


%% Face Centered
% FV solution average approximation (forward cell average at left face)
fprintf('FV approximation with a forward integral |     x-----|     |     |\n')
F = @(x) u(x);
Fh = @(x,dx) ubar(x,x+dx);
TE1 = @(x,dx) ux(x)/2*dx;
TE2 = @(x,dx) uxx(x)/6*dx^2;
TE3 = @(x,dx) uxxx(x)/24*dx^3;
TE4 = @(x,dx) uxxxx(x)/120*dx^4;
TE5 = @(x,dx) uxxxxx(x)/720*dx^5;
TE6 = @(x,dx) uxxxxxx(x)/5040*dx^6;

pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('FV solution average O(dx^1): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, TE1 , 10); fprintf('FV solution average O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) , 10); fprintf('FV solution average O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) , 9); fprintf('FV solution average O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) , 7); fprintf('FV solution average O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx) , 6); fprintf('FV solution average O(dx^6): %4.2f\n\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx) +TE6(x,dx) , 6); fprintf('FV solution average O(dx^7): %4.2f\n\n',pout);


% FV solution average approximation (forward cell average at left cell center)
fprintf('FV approximation with a forward integral at left cell center |     x     |-----| \n')
F = @(x) u(x);
Fh = @(x,dx) ubar(x+dx,x+2*dx);
TE1 = @(x,dx) ux(x)*3/2*dx;
TE2 = @(x,dx) uxx(x)*7/6*dx^2;
TE3 = @(x,dx) uxxx(x)*5/8*dx^3;
TE4 = @(x,dx) uxxxx(x)*31/120*dx^4;
TE5 = @(x,dx) uxxxxx(x)*7/80*dx^5;
TE6 = @(x,dx) uxxxxxx(x)*127/5040*dx^6;

pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('FV solution average O(dx^1): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, TE1 , 10); fprintf('FV solution average O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) , 10); fprintf('FV solution average O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) , 10); fprintf('FV solution average O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) , 8); fprintf('FV solution average O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx) , 6); fprintf('FV solution average O(dx^6): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx) +TE6(x,dx) , 6); fprintf('FV solution average O(dx^7): %4.2f\n\n',pout);


fprintf('FV approximation with a forward integral at left cell center |     x     |     |-----| \n')
F = @(x) u(x);
Fh = @(x,dx) ubar(x+2*dx,x+3*dx);
TE1 = @(x,dx) ux(x)*5/2*dx;
TE2 = @(x,dx) uxx(x)*19/6*dx^2;
TE3 = @(x,dx) uxxx(x)*65/24*dx^3;
TE4 = @(x,dx) uxxxx(x)*211/120*dx^4;
TE5 = @(x,dx) uxxxxx(x)*133/144*dx^5;
TE6 = @(x,dx) uxxxxxx(x)*2059/5040*dx^6;

pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('FV solution average O(dx^1): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, TE1 , 10); fprintf('FV solution average O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) , 10); fprintf('FV solution average O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) , 10); fprintf('FV solution average O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) , 8); fprintf('FV solution average O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx) , 6); fprintf('FV solution average O(dx^6): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx) +TE6(x,dx) , 6); fprintf('FV solution average O(dx^7): %4.2f\n\n',pout);


% FV solution average approximation (forward cell average at left cell center)
fprintf('FV approximation with a forward integral at left cell center  |     |  x  |-----| \n')
F = @(x) u(x);
Fh = @(x,dx) ubar(x+dx*1/2, x+dx*3/2);
TE1 = @(x,dx) ux(x)*dx;
TE2 = @(x,dx) uxx(x)*13/24*dx^2;
TE3 = @(x,dx) uxxx(x)*5/24*dx^3;
TE4 = @(x,dx) uxxxx(x)*121/1920*dx^4;
TE5 = @(x,dx) uxxxxx(x)*91/5760*dx^5;
TE6 = @(x,dx) uxxxxxx(x)*1093/322560*dx^6;

pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('FV solution average O(dx^1): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, TE1 , 10); fprintf('FV solution average O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) , 10); fprintf('FV solution average O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) , 10); fprintf('FV solution average O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) , 8); fprintf('FV solution average O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx) , 6); fprintf('FV solution average O(dx^6): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx) +TE6(x,dx) , 6); fprintf('FV solution average O(dx^7): %4.2f\n\n',pout);



fprintf('FV approximation with a forward integral at left cell center  |-----|  x  |     | \n')
F = @(x) u(x);
Fh = @(x,dx) ubar(x-dx*3/2, x-dx*1/2);
TE1 = @(x,dx) -ux(x)*dx;
TE2 = @(x,dx) uxx(x)*13/24*dx^2;
TE3 = @(x,dx) -uxxx(x)*5/24*dx^3;
TE4 = @(x,dx) uxxxx(x)*121/1920*dx^4;
TE5 = @(x,dx) -uxxxxx(x)*91/5760*dx^5;
TE6 = @(x,dx) uxxxxxx(x)*1093/322560*dx^6;

pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('FV solution average O(dx^1): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, TE1 , 10); fprintf('FV solution average O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) , 10); fprintf('FV solution average O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) , 10); fprintf('FV solution average O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) , 8); fprintf('FV solution average O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx) , 6); fprintf('FV solution average O(dx^6): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx) +TE6(x,dx) , 6); fprintf('FV solution average O(dx^7): %4.2f\n\n',pout);



fprintf('FV approximation with a forward integral at left cell center |     |  x  |     |-----| \n')
F = @(x) u(x);
Fh = @(x,dx) ubar(x+3*dx/2, x+dx*5/2);
TE1 = @(x,dx) ux(x)*2*dx;
TE2 = @(x,dx) uxx(x)*49/24*dx^2;
TE3 = @(x,dx) uxxx(x)*17/12*dx^3;
TE4 = @(x,dx) uxxxx(x)*1441/1920*dx^4;
TE5 = @(x,dx) uxxxxx(x)*931/2880*dx^5;
TE6 = @(x,dx) uxxxxxx(x)*37969/322560*dx^6;

pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('FV solution average O(dx^1): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, TE1 , 10); fprintf('FV solution average O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) , 10); fprintf('FV solution average O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) , 10); fprintf('FV solution average O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) , 8); fprintf('FV solution average O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx) , 6); fprintf('FV solution average O(dx^6): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx)+ TE6(x,dx) , 6); fprintf('FV solution average O(dx^7): %4.2f\n\n',pout);




% F = @(x) u(x).*ux(x);
% Fh = @(x,dx) (u(x+dx).^2 - u(x-dx).^2)/(4*dx);
% TE2 = @(x,dx) dx.^2*( 1/8*u(x)*uxxx(x) + 3/8*ux(x).*uxx(x));
% TE4 = @(x,dx) dx.^4*( u(x).*uxxxxx(x)/192 + 5/192*ux(x).*uxxxx(x) + 5/196*uxx(x)*uxxx(x));
% TE = @(x,dx) TE2(x,dx);
pout=OrderTest(x0, dx, F, Fh, TE, 10)



end

function pout=OrderTest(x, dx, F, Fh, TE, n)

    for i=1:n
       err(i) = abs(F(x) - Fh(x,dx) + TE(x,dx));
       TEex(i) = abs(F(x) - Fh(x,dx));
       TEa(i) = abs(TE(x,dx));

       dx = dx/2;


    end
    
    h = 1./2.^(0:n-1)./(1/2.^(n-1));
    p=log(err(1:end-1)./err(2:end))/log(2);

    loglog(h,err,'k-*',h,TEex,'r-o',h,TEa,'g-*')

    pout=p(end);
  
  
end

function [pout, p, err] = DerivativeTest(u,ux,x,dx,n)


    for i = 1:n
       f = ux(x);
       fh = (u(x+dx)-u(x-dx))/(2*dx);
       err(i) = fh-f;
       dx=dx/2;

    end

  p=log(err(1:end-1)./err(2:end))/log(2);
  
  pout=p(end);


end