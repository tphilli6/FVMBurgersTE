%% BTE_Neumann_inconsistent_4th


fprintf(' 4th Order extrapolation to face (centered at face):\n')

dudxA = @(x) u(x);
dudx = @(x,dx) -( 18*ux(x)*dx - 85*ubar(x,x+dx) + 23*ubar(x+dx,x+dx*2) - 4*ubar(x+dx*2,x+dx*3) )/66;
TE4 = @(x,dx) dx^4*(uxxxx(x)*3/110);
TE5 = @(x,dx) dx^5*(uxxxxx(x)*3/110);
TE6 = @(x,dx) dx^6*(uxxxxxx(x)*5/308);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) 0 , 10); fprintf('Extrapolate to face O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE4(x,dx) , 8); fprintf('Extrapolate to face O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE4(x,dx) + TE5(x,dx) , 8); fprintf('Extrapolate to face O(dx^6): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE4(x,dx) + TE5(x,dx) + TE6(x,dx), 7); fprintf('Extrapolate to face O(dx^7): %4.2f\n\n',pout);


fprintf(' 2nd Order BC flux at the face (centered at the face):\n')

F = @(x) u(x)^2/2 - nu*ux(x);
u0 = @(x,dx)  -( 18*ux(x)*dx - 85*ubar(x,x+dx) + 23*ubar(x+dx,x+dx*2) - 4*ubar(x+dx*2,x+dx*3) )/66;
Fh = @(x,dx) u0(x,dx).^2/2 - nu*ux(x);

TE4 = @(x,dx) dx^4*(u(x)*uxxxx(x)*3/110 + ux(x)*uxxx(x)*0 + uxx(x)^2*0/72  );
TE5 = @(x,dx) dx^5*(u(x)*uxxxxx(x)*3/110 + ux(x)*uxxxx(x)*0 + uxx(x)*uxxx(x)*0/144 );
TE6 = @(x,dx) dx^6*( u(x)*uxxxxxx(x)*5/308 );


pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('BC flux O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE4(x,dx), 7); fprintf('BC flux O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE4(x,dx) + TE5(x,dx), 7); fprintf('BC flux O(dx^6): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE4(x,dx) + TE5(x,dx) + TE6(x,dx), 6); fprintf('BC flux O(dx^7): %4.2f\n\n',pout);



fprintf(' 4th Order extrapolation to face (centered at cell center):\n')
TEcorr = @(x,dx) -1/2*ux(x)*dx ...
+ uxx(x)*1/8*dx^2 ...
- uxxx(x)*1/48*dx^3 ...
+ uxxxx(x)*1/384*dx^4 ...
- uxxxxx(x)*1/3840*dx^5 ...
+ uxxxxxx(x)*1/46080*dx^6;

dudxA = @(x) u(x);
dudx = @(x,dx) -( 18*ux(x-dx/2)*dx - 85*ubar(x-dx/2,x+dx/2) + 23*ubar(x+dx/2,x+dx*3/2) - 4*ubar(x+dx*3/2,x+dx*5/2) )/66 - TEcorr(x,dx);
TE4 = @(x,dx) dx^4*(uxxxx(x)*3/110);
TE5 = @(x,dx) dx^5*(uxxxxx(x)*3/220);
TE6 = @(x,dx) dx^6*(uxxxxxx(x)*37/6160);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) 0 , 10); fprintf('Extrapolate to face O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx)  TE4(x,dx), 7); fprintf('Extrapolate to face O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx)  TE4(x,dx) + TE5(x,dx), 6); fprintf('Extrapolate to face O(dx^6): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx)  TE4(x,dx) + TE5(x,dx) + TE6(x,dx), 6); fprintf('Extrapolate to face O(dx^7): %4.2f\n\n',pout);




fprintf(' 2nd Order BC flux at the face (centered at the cell center):\n')

F = @(x) u(x)^2/2 - nu*ux(x);
u0 = @(x,dx)  -( 18*ux(x-dx/2)*dx - 85*ubar(x-dx/2,x+dx/2) + 23*ubar(x+dx/2,x+dx*3/2) - 4*ubar(x+dx*3/2,x+dx*5/2) )/66;
Fh = @(x,dx) u0(x,dx).^2/2 - nu*ux(x-dx/2);

TE1 = @(x,dx) -dx*(u(x)*ux(x)*1/2 - uxx(x)*nu*1/2);
TE2 = @(x,dx) dx^2*( u(x)*uxx(x)*1/8 + ux(x)^2*1/8 - uxxx(x)*nu*1/8);
TE3 = @(x,dx) -dx^3*(u(x)*uxxx(x)*1/48 + ux(x)*uxx(x)*1/16 - uxxxx(x)*nu*1/48 );
TE4 = @(x,dx) dx^4*(u(x)*uxxxx(x)*631/21120 + ux(x)*uxxx(x)*1/96 + uxx(x)^2*1/128 - uxxxxx(x)*nu*1/384 );
TE5 = @(x,dx) -dx^5*(-u(x)*uxxxxx(x)*113/8448 + ux(x)*uxxxx(x)*631/42240 + uxx(x)*uxxx(x)*1/384 - uxxxxxx(x)*nu*1/3840 );


pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('BC flux O(dx^1): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx), 7); fprintf('BC flux O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx), 7); fprintf('BC flux O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx), 7); fprintf('BC flux O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx), 7); fprintf('BC flux O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx), 6); fprintf('BC flux O(dx^6): %4.2f\n\n',pout);



fprintf('Residual for boundary cell (centered at the cell center):\n')
F = @(x) u(x).*ux(x) - nu*uxx(x);
FhR = @(x,dx) ( ubar(x-dx*1/2,x+1/2*dx).^2+ubar(x+dx*1/2,x+dx*3/2).^2)/4 - nu*(ubar(x+dx*1/2,x+dx*3/2)-ubar(x-dx*1/2,x+1/2*dx))/dx;

u0 = @(x,dx)  -( 18*ux(x-dx/2)*dx - 85*ubar(x-dx/2,x+dx/2) + 23*ubar(x+dx/2,x+dx*3/2) - 4*ubar(x+dx*3/2,x+dx*5/2) )/66;
FhL = @(x,dx) u0(x,dx).^2/2- nu*ux(x-dx/2);

Fh = @(x,dx) (FhR(x,dx)-FhL(x,dx))/dx;

TE1 = @(x,dx) dx*( u(x)*uxx(x)*1/6 + ux(x)^2*1/8 - nu/12*uxxx(x));
TE2 = @(x,dx) dx.^2*( u(x).*uxxx(x)*1/8 + ux(x).*uxx(x)*1/3 - nu*uxxxx(x)*1/12); % + O(dx^4)
TE3 = @(x,dx) dx.^3*( u(x)*uxxxx(x)*1/528 + ux(x)*uxxx(x)*3/32 + uxx(x)^2*19/288 - nu*uxxxxx(x)*19/1440 );
TE4 = @(x,dx) dx.^4*( ux(x)*uxxxx(x)*327/7040 + uxx(x)*uxxx(x)*17/288 + u(x)*uxxxxx(x)*347/63360 - nu*uxxxxxx(x)*7/1920 ); % + O(dx^6)


pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('BC residual O(dx^1): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) , 9); fprintf('BC residual O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) , 8); fprintf('BC residual O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx), 7); fprintf('BC residual O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx), 7); fprintf('BC residual O(dx^5): %4.2f\n\n',pout);

