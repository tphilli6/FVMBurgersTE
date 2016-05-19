%% BTE_Dirichlet_inconsistent_3rd


fprintf(' 2nd Order extrapolation to face (centered at face):\n')
dudxA = @(x) ux(x);
dudx = @(x,dx) -( 6*u(x) - 7*ubar(x,x+dx) + ubar(x+dx,x+2*dx) )/(2*dx);
TE2 = @(x,dx) dx^2*(-uxxx(x)/6);
TE3 = @(x,dx) dx^3*(-uxxxx(x)/10);
TE4 = @(x,dx) dx^4*(-7/180*uxxxxx(x));
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) 0 , 10); fprintf('Extrapolate to face O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE2(x,dx) , 9); fprintf('Extrapolate to face O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE2(x,dx) + TE3(x,dx) , 8); fprintf('Extrapolate to face O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE2(x,dx) + TE3(x,dx) + TE4(x,dx), 7); fprintf('Extrapolate to face O(dx^5): %4.2f\n\n',pout);

fprintf(' 2nd Order BC flux at the face (centered at the face):\n')
F = @(x) u(x)^2/2 - nu*ux(x);
dudx = @(x,dx) -( 6*u(x) - 7*ubar(x,x+dx) + ubar(x+dx,x+2*dx) )/(2*dx);
Fh = @(x,dx) u(x)^2/2 - nu*dudx(x,dx);

TE2 = @(x,dx) dx^2*(nu*uxxx(x)*1/6);
TE3 = @(x,dx) dx^3*(nu*uxxxx(x)*1/10 );
TE4 = @(x,dx) dx^4*(nu*uxxxxx(x)*7/180 );
TE5 = @(x,dx) dx^5*(nu*uxxxxxx(x)*1/84 );


pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('BC flux O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) , 8); fprintf('BC flux O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) + TE3(x,dx), 7); fprintf('BC flux O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) + TE3(x,dx) + TE4(x,dx), 7); fprintf('BC flux O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx), 6); fprintf('BC flux O(dx^5): %4.2f\n\n',pout);



fprintf(' 2nd Order extrapolation to face (centered at cell center):\n')
dudxA = @(x) -ux(x);
dudx = @(x,dx) -(-( 6*u(x-dx/2) - 7*ubar(x-dx/2,x+dx/2) + ubar(x+dx/2,x+3/2*dx) )/(2*dx));
TE1 = @(x,dx) dx*(uxx(x)/2);
TE2 = @(x,dx) dx^2*(uxxx(x)/24);
TE3 = @(x,dx) dx^3*(uxxxx(x)*3/80);
TE4 = @(x,dx) dx^4*(uxxxxx(x)*41/5760);
TE5 = @(x,dx) dx^5*(uxxxxxx(x)*47/26880);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) 0 , 10); fprintf('Extrapolate to face O(dx^1): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE1(x,dx) , 9); fprintf('Extrapolate to face O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE1(x,dx) + TE2(x,dx) , 8); fprintf('Extrapolate to face O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx), 7); fprintf('Extrapolate to face O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx), 6); fprintf('Extrapolate to face O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx), 5); fprintf('Extrapolate to face O(dx^6): %4.2f\n\n',pout);




fprintf(' 2nd Order BC flux at the face (centered at the cell center):\n')
F = @(x) u(x)^2/2 - nu*ux(x);
dudx = @(x,dx) -(6*u(x-dx/2) - 7*ubar(x-dx/2,x+dx/2) + ubar(x+dx/2,x+3/2*dx))/(2*dx);
Fh = @(x,dx) u(x-dx/2)^2/2 - nu*dudx(x,dx);

TE1 = @(x,dx) -dx*(u(x)*ux(x)/2 - nu*uxx(x)/2 );
TE2 = @(x,dx) dx^2*(ux(x)^2/8 + u(x)*uxx(x)/8 + nu*uxxx(x)/24);
TE3 = @(x,dx) -dx^3*(u(x)*uxxx(x)*1/48 + ux(x)*uxx(x)/16 - nu*uxxxx(x)*3/80 );
TE4 = @(x,dx) dx^4*(u(x)*uxxxx(x)/384 + ux(x)*uxxx(x)/96 + uxx(x)^2/128 + nu*uxxxxx(x)*41/5760 );
TE5 = @(x,dx) -dx^5*(u(x)*uxxxxx(x)*1/3840 + ux(x)*uxxxx(x)/768 + uxx(x)*uxxx(x)*1/384 - nu*uxxxxxx(x)*47/26880 );


pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('BC flux O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) , 9); fprintf('BC flux O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) , 8); fprintf('BC flux O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx), 7); fprintf('BC flux O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx), 7); fprintf('BC flux O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx), 6); fprintf('BC flux O(dx^5): %4.2f\n\n',pout);



fprintf('Residual for boundary cell (centered at the cell center):\n')
F = @(x) u(x).*ux(x) - nu*uxx(x);
FhR = @(x,dx) ( ubar(x-dx*1/2,x+1/2*dx).^2+ubar(x+dx*1/2,x+dx*3/2).^2)/4 - nu*(ubar(x+dx*1/2,x+dx*3/2)-ubar(x-dx*1/2,x+1/2*dx))/dx;

dudx = @(x,dx) -(6*u(x-dx/2) - 7*ubar(x-dx/2,x+dx/2) + ubar(x+dx/2,x+3/2*dx))/(2*dx);
FhL = @(x,dx) u(x-dx/2)^2/2 - nu*dudx(x,dx);
Fh = @(x,dx) (FhR(x,dx)-FhL(x,dx))/dx;

TE1 = @(x,dx) dx*( 1/8*ux(x)^2 + 1/6*u(x)*uxx(x) - nu/4*uxxx(x));
TE2 = @(x,dx) dx.^2*( u(x).*uxxx(x)*1/8 + ux(x).*uxx(x)*1/3 - nu*uxxxx(x)/10); % + O(dx^4)
TE3 = @(x,dx) dx.^3*( u(x)*uxxxx(x)*7/240 + uxx(x)^2*19/288 + ux(x)*uxxx(x)*3/32 - nu*uxxxxx(x)*11/480 );
TE4 = @(x,dx) dx.^4*( ux(x)*uxxxx(x)*21/640 + uxx(x)*uxxx(x)*17/288 + u(x)*uxxxxx(x)*47/5760 - nu*uxxxxxx(x)*23/4480 ); % + O(dx^6)


pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('BC residual O(dx^1): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) , 9); fprintf('BC residual O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) , 8); fprintf('BC residual O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx), 7); fprintf('BC residual O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx), 6); fprintf('BC residual O(dx^5): %4.2f\n\n',pout);

