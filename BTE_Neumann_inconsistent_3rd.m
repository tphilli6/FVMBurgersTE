%% BTE_Neumann_inconsistent_3rd


fprintf(' 3rd Order extrapolation to face (centered at face):\n')

dudxA = @(x) u(x);
dudx = @(x,dx) -( 2*ux(x)*dx - 7*ubar(x,x+dx) + ubar(x+dx,x+dx*2))/6;
TE3 = @(x,dx) -dx^3*(uxxx(x)/18);
TE4 = @(x,dx) -dx^4*(uxxxx(x)*1/30);
TE5 = @(x,dx) -dx^5*(7/540*uxxxxx(x));
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) 0 , 10); fprintf('Extrapolate to face O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE3(x,dx) , 8); fprintf('Extrapolate to face O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE3(x,dx) + TE4(x,dx) , 8); fprintf('Extrapolate to face O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) TE3(x,dx) + TE4(x,dx) + TE5(x,dx), 7); fprintf('Extrapolate to face O(dx^6): %4.2f\n\n',pout);


fprintf(' 2nd Order BC flux at the face (centered at the face):\n')

F = @(x) u(x)^2/2 - nu*ux(x);
u0 = @(x,dx) -( 2*ux(x)*dx - 7*ubar(x,x+dx) + ubar(x+dx,x+dx*2))/6;
Fh = @(x,dx) u0(x,dx).^2/2 - nu*ux(x);

TE2 = @(x,dx) 0;%dx^2*(ux(x)^2*1/8 + u(x)*uxx(x)*1/6);
TE3 = @(x,dx) dx^3*(-u(x)*uxxx(x)*1/18  );
TE4 = @(x,dx) dx^4*(-u(x)*uxxxx(x)*1/30 + ux(x)*uxxx(x)*0 + uxx(x)^2*0/72  );
TE5 = @(x,dx) dx^5*(-u(x)*uxxxxx(x)*7/540 + ux(x)*uxxxx(x)*0 + uxx(x)*uxxx(x)*0/144 );



pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('BC flux O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) + TE3(x,dx), 7); fprintf('BC flux O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) + TE3(x,dx) + TE4(x,dx), 7); fprintf('BC flux O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE2(x,dx) + TE3(x,dx) + TE4(x,dx) + TE5(x,dx), 6); fprintf('BC flux O(dx^6): %4.2f\n\n',pout);



fprintf(' 2nd Order extrapolation to face (centered at cell center):\n')
TEcorr = @(x,dx) -1/2*ux(x)*dx ...
+ uxx(x)*1/8*dx^2 ...
- uxxx(x)*1/48*dx^3 ...
+ uxxxx(x)*1/384*dx^4 ...
- uxxxxx(x)*1/3840*dx^5 ...
+ uxxxxxx(x)*1/46080*dx^6;

dudxA = @(x) u(x);

u0 = @(x,dx) -( 2*ux(x)*dx - 7*ubar(x,x+dx) + ubar(x+dx,x+dx*2))/6;
dudx = @(x,dx) -( 2*ux(x-dx/2)*dx - 7*ubar(x-dx/2,x+dx/2) + ubar(x+dx/2,x+dx*3/2))/6 - TEcorr(x,dx);
TE3 = @(x,dx) -dx^3*(uxxx(x)*1/18);
TE4 = @(x,dx) -dx^4*(uxxxx(x)*1/180);
TE5 = @(x,dx) -dx^5*(uxxxxx(x)*7/2160);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx) 0 , 10); fprintf('Extrapolate to face O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx)  TE3(x,dx), 7); fprintf('Extrapolate to face O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx)  TE3(x,dx) + TE4(x,dx), 6); fprintf('Extrapolate to face O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, dudxA, dudx, @(x,dx)  TE3(x,dx) + TE4(x,dx) + TE5(x,dx), 6); fprintf('Extrapolate to face O(dx^6): %4.2f\n\n',pout);




fprintf(' 2nd Order BC flux at the face (centered at the cell center):\n')

F = @(x) u(x)^2/2 - nu*ux(x);
u0 = @(x,dx)  -( 2*ux(x-dx/2)*dx - 7*ubar(x-dx/2,x+dx/2) + ubar(x+dx/2,x+dx*3/2))/6;
Fh = @(x,dx) u0(x,dx).^2/2- nu*ux(x-dx/2);

TE3 = @(x,dx) -dx^3*(u(x)*uxxx(x)*11/144 + ux(x)*uxx(x)*1/16 - nu*uxxxx(x)*1/48 );
TE4 = @(x,dx) dx^4*(-u(x)*uxxxx(x)*17/5760 + ux(x)*uxxx(x)*11/288 + uxx(x)^2*1/128 - nu*uxxxxx(x)*1/384 );
TE5 = @(x,dx) -dx^5*(u(x)*uxxxxx(x)*121/34560 - ux(x)*uxxxx(x)*17/11520 + uxx(x)*uxxx(x)*11/1152 - nu*uxxxxxx(x)*1/3840 );


pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('BC flux O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE3(x,dx), 7); fprintf('BC flux O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE3(x,dx) + TE4(x,dx), 7); fprintf('BC flux O(dx^5): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE3(x,dx) + TE4(x,dx) + TE5(x,dx), 6); fprintf('BC flux O(dx^6): %4.2f\n\n',pout);



fprintf('Residual for boundary cell (centered at the cell center):\n')
F = @(x) u(x).*ux(x) - nu*uxx(x);
FhR = @(x,dx) ( ubar(x-dx*1/2,x+1/2*dx).^2+ubar(x+dx*1/2,x+dx*3/2).^2)/4 - nu*(ubar(x+dx*1/2,x+dx*3/2)-ubar(x-dx*1/2,x+1/2*dx))/dx;

u0 = @(x,dx)  -( 2*ux(x-dx/2)*dx - 7*ubar(x-dx/2,x+dx/2) + ubar(x+dx/2,x+dx*3/2))/6;
FhL = @(x,dx) u0(x,dx).^2/2- nu*ux(x-dx/2);

Fh = @(x,dx) (FhR(x,dx)-FhL(x,dx))/dx;

TE1 = @(x,dx) dx*( u(x)*uxx(x)*1/6 + ux(x)^2*1/8 - nu/12*uxxx(x));
TE2 = @(x,dx) dx.^2*( u(x).*uxxx(x)*13/72 + ux(x).*uxx(x)*1/3 - nu*uxxxx(x)*1/12); % + O(dx^4)
TE3 = @(x,dx) dx.^3*( u(x)*uxxxx(x)*5/144 + ux(x)*uxxx(x)*19/288 + uxx(x)^2*19/288 - nu*uxxxxx(x)*19/1440 );
TE4 = @(x,dx) dx.^4*( ux(x)*uxxxx(x)*173/5760 + uxx(x)*uxxx(x)*19/288 + u(x)*uxxxxx(x)*197/17280 - nu*uxxxxxx(x)*7/1920 ); % + O(dx^6)


pout=OrderTest(x0, dx, F, Fh, @(x,dx) 0 , 10); fprintf('BC residual O(dx^1): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) , 9); fprintf('BC residual O(dx^2): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) , 8); fprintf('BC residual O(dx^3): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx), 7); fprintf('BC residual O(dx^4): %4.2f\n',pout);
pout=OrderTest(x0, dx, F, Fh, @(x,dx) TE1(x,dx) + TE2(x,dx) + TE3(x,dx) + TE4(x,dx), 7); fprintf('BC residual O(dx^5): %4.2f\n\n',pout);

