%% BTE_Neumann_consistent_3rd


fprintf(' 3rd Order extrapolation to ghost cell (centered at face):\n')
TEcorr = @(x,dx) -ux(x)/2*dx ...
 +uxx(x)/6*dx^2 ...
 -uxxx(x)/24*dx^3 ...
 +uxxxx(x)/120*dx^4 ...
 -uxxxxx(x)/720*dx^5;

dudxA = @(x) u(x);
dudx = @(x,dx) ubar(x,x+dx) - ux(x)*dx - TEcorr(x,dx);
TE3 = @(x,dx) dx^3*(uxxx(x)/12);
TE5 = @(x,dx) dx^5*(1/360*uxxxxx(x));
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) 0 , 10); fprintf('Extrapolate to face O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE3(x,dx) , 8); fprintf('Extrapolate to face O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE3(x,dx)+ TE5(x,dx), 5); fprintf('Extrapolate to face O(dx^7): %4.2f\n\n',pout);


fprintf(' 2nd Order BC flux at the face (centered at the face):\n')

F = @(x) u(x)^2/2 - nu*ux(x);
u0 = @(x,dx) ubar(x,x+dx) - ux(x)*dx;
Fh = @(x,dx) (ubar(x,x+dx)^2 + u0(x,dx).^2)/4 - nu*(ubar(x,x+dx)-u0(x,dx))/dx;

TE2 = @(x,dx) dx^2*(ux(x)^2*1/8 + u(x)*uxx(x)*1/6);
TE3 = @(x,dx) dx^3*(u(x)*uxxx(x)*1/24  );
TE4 = @(x,dx) dx^4*(u(x)*uxxxx(x)*1/120 + ux(x)*uxxx(x)*0 + uxx(x)^2*1/72  );
TE5 = @(x,dx) dx^5*(u(x)*uxxxxx(x)*1/720 + ux(x)*uxxxx(x)*0 + uxx(x)*uxxx(x)*1/144 );



pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('BC flux O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) , 8); fprintf('BC flux O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) + TE3(x,dx), 7); fprintf('BC flux O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) + TE3(x,dx) + TE4(x,dx), 7); fprintf('BC flux O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx), 6); fprintf('BC flux O(dx^6): %4.2f\n\n',pout);



fprintf(' 2nd Order extrapolation to face (centered at cell center):\n')
TEcorr = @(x,dx) -ux(x)*dx ...
+ uxx(x)*13/24*dx^2 ...
- uxxx(x)*5/24*dx^3 ...
+ uxxxx(x)*121/1920*dx^4 ...
- uxxxxx(x)*91/5760*dx^5 ...
+ uxxxxxx(x)*1093/322560*dx^6;

dudxA = @(x) u(x);
dudx = @(x,dx) ubar(x-dx/2,x+dx/2) - ux(x-dx/2)*dx - TEcorr(x,dx);
TE3 = @(x,dx) dx^3*(uxxx(x)*1/12);
TE4 = @(x,dx) -dx^4*(uxxxx(x)*1/24);
TE5 = @(x,dx) dx^5*(uxxxxx(x)*19/1440);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) 0 , 10); fprintf('Extrapolate to face O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx)  TE3(x,dx), 7); fprintf('Extrapolate to face O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx)  TE3(x,dx) + TE4(x,dx), 6); fprintf('Extrapolate to face O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx)  TE3(x,dx) + TE4(x,dx) + TE5(x,dx), 5); fprintf('Extrapolate to face O(dx^6): %4.2f\n\n',pout);




fprintf(' 2nd Order BC flux at the face (centered at the cell center):\n')

F = @(x) u(x)^2/2 - nu*ux(x);
u0 = @(x,dx) ubar(x-dx/2,x+dx/2) - ux(x-dx/2)*dx;
Fh = @(x,dx) (ubar(x-dx/2,x+dx/2)^2 + u0(x,dx).^2)/4 - nu*ux(x-dx/2);

TE1 = @(x,dx) -dx*(u(x)*ux(x)/2 - nu*uxx(x)/2 );
TE2 = @(x,dx) dx^2*(ux(x)^2/4 + u(x)*uxx(x)*7/24 - nu*uxxx(x)/8);
TE3 = @(x,dx) -dx^3*(u(x)*uxxx(x)*1/16 + ux(x)*uxx(x)*13/48 - nu*uxxxx(x)*1/48 );
TE4 = @(x,dx) dx^4*(u(x)*uxxxx(x)*7/640 + ux(x)*uxxx(x)*1/16 + uxx(x)^2*85/1152 + nu*uxxxxx(x)*1/384 );
TE5 = @(x,dx) -dx^5*(u(x)*uxxxxx(x)*1/768 + ux(x)*uxxxx(x)*41/3840 + uxx(x)*uxxx(x)*13/384 - nu*uxxxxxx(x)*1/3840 );


pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('BC flux O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) , 9); fprintf('BC flux O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) , 8); fprintf('BC flux O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx), 7); fprintf('BC flux O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx), 7); fprintf('BC flux O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx), 6); fprintf('BC flux O(dx^6): %4.2f\n\n',pout);



fprintf('Residual for boundary cell (centered at the cell center):\n')
F = @(x) u(x).*ux(x) - nu*uxx(x);
FhR = @(x,dx) ( ubar(x-dx*1/2,x+1/2*dx).^2+ubar(x+dx*1/2,x+dx*3/2).^2)/4 - nu*(ubar(x+dx*1/2,x+dx*3/2)-ubar(x-dx*1/2,x+1/2*dx))/dx;

u0 = @(x,dx) ubar(x-dx/2,x+dx/2) - ux(x-dx/2)*dx;
FhL = @(x,dx) (ubar(x-dx/2,x+dx/2)^2 + u0(x,dx).^2)/4 - nu*(ubar(x-dx/2,x+dx/2)-u0(x,dx))/dx;

Fh = @(x,dx) (FhR(x,dx)-FhL(x,dx))/dx;

TE1 = @(x,dx) dx*( - nu/12*uxxx(x));
TE2 = @(x,dx) dx.^2*( u(x).*uxxx(x)*1/6 + ux(x).*uxx(x)*13/24 - nu*uxxxx(x)*1/12); % + O(dx^4)
TE3 = @(x,dx) dx.^3*( u(x)*uxxxx(x)*1/48 + ux(x)*uxxx(x)*1/24 - nu*uxxxxx(x)*19/1440 );
TE4 = @(x,dx) dx.^4*( ux(x)*uxxxx(x)*27/640 + uxx(x)*uxxx(x)*13/144 + u(x)*uxxxxx(x)*53/5760 - nu*uxxxxxx(x)*7/1920 ); % + O(dx^6)


pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('BC residual O(dx^1): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) , 9); fprintf('BC residual O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) , 8); fprintf('BC residual O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx), 7); fprintf('BC residual O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx), 7); fprintf('BC residual O(dx^5): %4.2f\n\n',pout);

