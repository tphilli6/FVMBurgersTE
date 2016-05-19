%% BTE_Neumann_consistent_4th


fprintf(' 3rd Order extrapolation to ghost cell (centered at face):\n')
TEcorr = @(x,dx) -ux(x)/2*dx ...
 +uxx(x)/6*dx^2 ...
 -uxxx(x)/24*dx^3 ...
 +uxxxx(x)/120*dx^4 ...
 -uxxxxx(x)/720*dx^5;

dudxA = @(x) u(x);
dudx = @(x,dx) (-12*ux(x)*dx + 9*ubar(x,x+dx) + 3*ubar(x+dx,x+2*dx) - ubar(x+dx*2,x+dx*3))/11- TEcorr(x,dx);
TE4 = @(x,dx) -dx^4*(uxxxx(x)/11);
TE5 = @(x,dx) -dx^5*(19/330*uxxxxx(x));
TE6 = @(x,dx) -dx^6*(1/33*uxxxxxx(x));
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) 0 , 10); fprintf('Extrapolate to face O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE4(x,dx) , 8); fprintf('Extrapolate to face O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE4(x,dx)+ TE5(x,dx), 5); fprintf('Extrapolate to face O(dx^6): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE4(x,dx)+ TE5(x,dx)+ TE6(x,dx), 6); fprintf('Extrapolate to face O(dx^7): %4.2f\n\n',pout);

fprintf(' 2nd Order BC flux at the face (centered at the face):\n')

F = @(x) u(x)^2/2 - nu*ux(x);
u0 = @(x,dx)  (-12*ux(x)*dx + 9*ubar(x,x+dx) + 3*ubar(x+dx,x+2*dx) - ubar(x+dx*2,x+dx*3))/11;
Fh = @(x,dx) (ubar(x,x+dx)^2 + u0(x,dx).^2)/4 - nu*(ubar(x,x+dx)-u0(x,dx))/dx;

TE2 = @(x,dx) dx^2*(ux(x)^2*1/8 + u(x)*uxx(x)*1/6 - nu*uxxx(x)*1/12);
TE3 = @(x,dx) -dx^3*(nu*uxxxx(x)*1/11 );
TE4 = @(x,dx) dx^4*(-u(x)*uxxxx(x)*49/1320 + ux(x)*uxxx(x)*1/48 + uxx(x)^2*1/72 - nu*uxxxxx(x)*239/3960 );
TE5 = @(x,dx) +dx^5*(-u(x)*uxxxxx(x)*19/660 + ux(x)*uxxxx(x)*1/44 - uxx(x)*uxxx(x)*0 - nu*uxxxxxx(x)*1/33 );


pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('BC flux O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx)  TE2(x,dx) , 8); fprintf('BC flux O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx)  TE2(x,dx) + TE3(x,dx), 7); fprintf('BC flux O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx)  TE2(x,dx) + TE3(x,dx) + TE4(x,dx), 7); fprintf('BC flux O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx)  TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx), 6); fprintf('BC flux O(dx^6): %4.2f\n\n',pout);



fprintf(' 2nd Order extrapolation to face (centered at cell center):\n')
TEcorr = @(x,dx) -ux(x)*dx ...
+ uxx(x)*13/24*dx^2 ...
- uxxx(x)*5/24*dx^3 ...
+ uxxxx(x)*121/1920*dx^4 ...
- uxxxxx(x)*91/5760*dx^5 ...
+ uxxxxxx(x)*1093/322560*dx^6;

dudxA = @(x) u(x);
dudx = @(x,dx) (-12*ux(x-dx/2)*dx + 9*ubar(x-dx/2,x+dx/2) + 3*ubar(x+dx/2,x+dx*3/2) - ubar(x+dx*3/2,x+dx*5/2))/11- TEcorr(x,dx);
TE4 = @(x,dx) -dx^4*(uxxxx(x)*1/11);
TE5 = @(x,dx) -dx^5*(uxxxxx(x)*2/165);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) 0 , 10); fprintf('Extrapolate to face O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx)   TE4(x,dx), 6); fprintf('Extrapolate to face O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx)   TE4(x,dx) + TE5(x,dx), 7); fprintf('Extrapolate to face O(dx^6): %4.2f\n\n',pout);




fprintf(' 2nd Order BC flux at the face (centered at the cell center):\n')

F = @(x) u(x)^2/2 - nu*ux(x);
u0 = @(x,dx) (-12*ux(x-dx/2)*dx + 9*ubar(x-dx/2,x+dx/2) + 3*ubar(x+dx/2,x+dx*3/2) - ubar(x+dx*3/2,x+dx*5/2))/11;

Fh = @(x,dx) (ubar(x-dx/2,x+dx/2)^2 + u0(x,dx).^2)/4 - nu*(ubar(x-dx/2,x+dx/2)-u0(x,dx))/dx;

TE2 = @(x,dx) dx^2*(ux(x)^2*1/4 + u(x)*uxx(x)*7/24 - nu*uxxx(x)*5/24);
TE3 = @(x,dx) -dx^3*(u(x)*uxxx(x)*5/48 + ux(x)*uxx(x)*13/48 + nu*uxxxx(x)*5/176 );
TE4 = @(x,dx) dx^4*(-u(x)*uxxxx(x)*289/21120 + ux(x)*uxxx(x)*5/48 + uxx(x)^2*85/1152 - nu*uxxxxx(x)*1769/63360 );
TE5 = @(x,dx) +dx^5*(-u(x)*uxxxxx(x)*1769/126720 + ux(x)*uxxxx(x)*589/42240 - uxx(x)*uxxx(x)*65/1152 - nu*uxxxxxx(x)*101/10639 );


pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('BC flux O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx)  TE2(x,dx) , 8); fprintf('BC flux O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx)  TE2(x,dx) + TE3(x,dx), 7); fprintf('BC flux O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx)  TE2(x,dx) + TE3(x,dx) + TE4(x,dx), 7); fprintf('BC flux O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx)  TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx), 6); fprintf('BC flux O(dx^6): %4.2f\n\n',pout);



fprintf('Residual for boundary cell (centered at the cell center):\n')
F = @(x) u(x).*ux(x) - nu*uxx(x);
FhR = @(x,dx) ( ubar(x-dx*1/2,x+1/2*dx).^2+ubar(x+dx*1/2,x+dx*3/2).^2)/4 - nu*(ubar(x+dx*1/2,x+dx*3/2)-ubar(x-dx*1/2,x+1/2*dx))/dx;

u0 = @(x,dx) (-12*ux(x-dx/2)*dx + 9*ubar(x-dx/2,x+dx/2) + 3*ubar(x+dx/2,x+dx*3/2) - ubar(x+dx*3/2,x+dx*5/2))/11;
FhL = @(x,dx) (ubar(x-dx/2,x+dx/2)^2 + u0(x,dx).^2)/4 - nu*(ubar(x-dx/2,x+dx/2)-u0(x,dx))/dx

Fh = @(x,dx) (FhR(x,dx)-FhL(x,dx))/dx;

TE2 = @(x,dx) dx.^2*( u(x).*uxxx(x)*5/24 + ux(x).*uxx(x)*13/24 - nu*uxxxx(x)*3/88); % + O(dx^4)
TE3 = @(x,dx) dx.^3*( u(x)*uxxxx(x)*1/22 + ux(x)*uxxx(x)*0 + nu*uxxxxx(x)*2/165 );
TE4 = @(x,dx) dx.^4*( ux(x)*uxxxx(x)*372/21120 + uxx(x)*uxxx(x)*65/576 + u(x)*uxxxxx(x)*277/12672 + nu*uxxxxxx(x)*42/7040 ); % + O(dx^6)


pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('BC residual O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) , 8); fprintf('BC residual O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) + TE3(x,dx), 7); fprintf('BC residual O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) + TE3(x,dx) + TE4(x,dx), 7); fprintf('BC residual O(dx^5): %4.2f\n\n',pout);

