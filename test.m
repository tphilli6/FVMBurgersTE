clear all

syms nu dx u ux uxx uxxx uxxxx uxxxxx uxxxxxx dx1 dx2 dx3 reals %ut uxt uxxt uxxxt uxxxxt uxxxxxt uxxxxxxt c1 c2 c3 f fx fxx fxxx fxxxx fxxxxx nu ubar reals
clc


uf_im1 = u - 3/2*ux*dx +  7/6*uxx*dx^2  -   5/8*uxxx*dx^3  +   31/120*uxxxx*dx^4  -    7/80*uxxxxx*dx^5  +  127/5040*uxxxxxx*dx^6;
uf_i   = u - 1/2*ux*dx +  1/6*uxx*dx^2  -  1/24*uxxx*dx^3  +    1/120*uxxxx*dx^4  -   1/720*uxxxxx*dx^5  +    1/5040*uxxxxxx*dx^6;
uf_ip1 = u + 1/2*ux*dx +  1/6*uxx*dx^2  +  1/24*uxxx*dx^3  +    1/120*uxxxx*dx^4  +   1/720*uxxxxx*dx^5  +    1/5040*uxxxxxx*dx^6;
uf_ip2 = u + 3/2*ux*dx +  7/6*uxx*dx^2  +   5/8*uxxx*dx^3  +   31/120*uxxxx*dx^4  +    7/80*uxxxxx*dx^5  +  127/5040*uxxxxxx*dx^6;
uf_ip3 = u + 5/2*ux*dx + 19/6*uxx*dx^2  + 65/24*uxxx*dx^3  +  211/120*uxxxx*dx^4  + 133/144*uxxxxx*dx^5  + 2059/5040*uxxxxxx*dx^6;

uc_im3 = u - 3*ux*dx + 109/24*uxx*dx^2 -  37/8*uxxx*dx^3 + 6841/1920*uxxxx*dx^4 - 1417/2880*uxxxxx*dx^5 + 372709/322560*uxxxxxx*dx^6;
uc_im2 = u - 2*ux*dx +  49/24*uxx*dx^2 - 17/12*uxxx*dx^3 + 1441/1920*uxxxx*dx^4 -  931/2880*uxxxxx*dx^5 +  37969/322560*uxxxxxx*dx^6;
uc_im1 = u -   ux*dx +  13/24*uxx*dx^2 -  5/24*uxxx*dx^3 +  121/1920*uxxxx*dx^4 -   91/5760*uxxxxx*dx^5 +   1093/322560*uxxxxxx*dx^6;
uc_i   = u           +   1/24*uxx*dx^2                   +    1/1920*uxxxx*dx^4                         +      1/322560*uxxxxxx*dx^6;
uc_ip1 = u +   ux*dx +  13/24*uxx*dx^2 +  5/24*uxxx*dx^3 +  121/1920*uxxxx*dx^4 +   91/5760*uxxxxx*dx^5 +   1093/322560*uxxxxxx*dx^6;
uc_ip2 = u + 2*ux*dx +  49/24*uxx*dx^2 + 17/12*uxxx*dx^3 + 1441/1920*uxxxx*dx^4 +  931/2880*uxxxxx*dx^5 +  37969/322560*uxxxxxx*dx^6;
uc_ip3 = u + 3*ux*dx + 109/24*uxx*dx^2 +  37/8*uxxx*dx^3 + 6841/1920*uxxxx*dx^4 + 1417/2880*uxxxxx*dx^5 + 372709/322560*uxxxxxx*dx^6;

uxc_iph =  ux + 1/2*uxx*dx  +  1/8*uxxx*dx^2 +   1/48*uxxxx*dx^3  +  1/384*uxxxxx*dx^4 +    1/3840*uxxxxxx*dx^5;
uxc_ip3h = ux + 3/2*uxx*dx  +  9/8*uxxx*dx^2 +  27/48*uxxxx*dx^3  + 81/384*uxxxxx*dx^4 + 2431/3840*uxxxxxx*dx^5;
uxc_ip5h = ux + 5/3*uxx*dx  + 25/8*uxxx*dx^2 + 125/48*uxxxx*dx^3  +625/384*uxxxxx*dx^4 + 3125/3840*uxxxxxx*dx^5;


u_p1 = u  + ux*dx +  1/2*uxx*dx^2  +  1/6*uxxx*dx^3  +    1/24*uxxxx*dx^4  +  1/factorial(5)*uxxxxx*dx^5;
u_p2 = u  + ux*dx*2 +  4/2*uxx*dx^2  +  8/6*uxxx*dx^3  +    16/24*uxxxx*dx^4  +  32/factorial(5)*uxxxxx*dx^5;
u_p3 = u  + ux*dx*3 +  9/2*uxx*dx^2  +  27/6*uxxx*dx^3  +    81/24*uxxxx*dx^4  +  243/factorial(5)*uxxxxx*dx^5;

u_m1 = u  - ux*dx   +  1/2*uxx*dx^2  -  1/6*uxxx*dx^3  +    1/24*uxxxx*dx^4  +  1/factorial(5)*uxxxxx*dx^5;
u_m2 = u  - ux*dx*2 +  4/2*uxx*dx^2  -  8/6*uxxx*dx^3  +    16/24*uxxxx*dx^4  +  32/factorial(5)*uxxxxx*dx^5;
u_m3 = u  - ux*dx*3 +  9/2*uxx*dx^2  -  27/6*uxxx*dx^3  +    81/24*uxxxx*dx^4  +  243/factorial(5)*uxxxxx*dx^5;

u_m1t = u  - ux*(dx1)   +  1/2*uxx*(dx1)^2  -  1/6*uxxx*(dx1)^3;%  +    1/24*uxxxx*(dx1)^4;%  +  1/factorial(5)*uxxxxx*(dx1)^5;
u_m2t = u  - ux*(dx1+dx2)   +  1/2*uxx*(dx1+dx2)^2  -  1/6*uxxx*(dx1+dx2)^3;% +     1/24*uxxxx*(dx1+dx2)^4;%  +  1/factorial(5)*uxxxxx*(dx1+dx2)^5;
u_m3t = u  - ux*(dx1+dx2+dx3)   +  1/2*uxx*(dx1+dx2+dx3)^2  -  1/6*uxxx*(dx1+dx2+dx3)^3;%  +    1/24*uxxxx*(dx1+dx2+dx3)^4;%  +  1/factorial(5)*uxxxxx*(dx1+dx2+dx3)^5;


u_mh = u - 1/2*ux*dx +  1/8*uxx*dx^2  -  1/48*uxxx*dx^3  +    1/384*uxxxx*dx^4  -   1/3840*uxxxxx*dx^5  +    1/46080*uxxxxxx*dx^6;
ux_mh = ux - 1/2*uxx*dx +  1/8*uxxx*dx^2  -  1/48*uxxxx*dx^3  +    1/384*uxxxxx*dx^4  -   1/3840*uxxxxxx*dx^5;



dfdx = (3*u - 4*u_m1 + u_m2)/(2*dx);
fnph = u + 1/2*(u - u_m1);
fnmh = u_m1 + 1/2*(u_m1 - u_m2);


pretty(expand(fnph))
pretty(expand(fnmh))
dfdx2 = (fnph-fnmh)/dx;
pretty(expand(dfdx2));

dfdx2_nonequal = (u-u_m1t)/dx1 + (u-u_m1t)/(dx1+dx2) - dx2/dx1 * (u_m1t - u_m2t)/(dx2+dx3);
pretty(expand(dfdx2_nonequal))
pretty(expand(subs(dfdx2_nonequal,{dx1, dx2, dx3}, {dx, dx, dx}) ));
% fnmh

%% Face Centered Flux (Fh)
fh_fiph = (uf_ip1^2 + uf_i^2)/4 - nu*(uf_ip1-uf_i)/dx;
fh_fimh = (uf_i^2 + uf_im1^2)/4 - nu*(uf_i-uf_im1)/dx;

pretty(expand(simplify(fh_fiph)))


fh_ciph = (uc_ip1^2 + uc_i^2)/4 - nu*(uc_ip1-uc_i)/dx;
fh_cimh = (uc_i^2 + uc_im1^2)/4 - nu*(uc_i-uc_im1)/dx;

pretty(expand(simplify(fh_ciph)))
pretty(expand(simplify(fh_cimh)))
pretty(expand( (fh_ciph-fh_cimh)/dx ))


%% Scheme inconsistent Dirichlet
    %2nd order extrap 
    uxf = -(6*u - 7*uf_ip1 + uf_ip2)/(2*dx); %Face centered
%     pretty(simplify(expand(uxf)))

    ux_h = u^2/2 - nu*(uxf); %Face centered
%     pretty(simplify(expand(ux_h)))

    fh_cbc =u_mh^2/2 - nu*( -( 6*u_mh - 7*uc_i + uc_ip1)/(2*dx));%cell centered
%     pretty(simplify(expand(fh_cbc)))

    R_i1 = (fh_ciph-fh_cbc)/dx;
%     pretty(simplify(expand( R_i1 )))

% 4th order (3rd order derivative
    uxf = -(66*u - 85*uf_ip1 + 23*uf_ip2 - 4*uf_ip3)/(18*dx);
    pretty(expand(uxf))
    
    ux_h = u^2/2 - nu*uxf; %Face centered
    pretty(simplify(expand(ux_h)))
   
    uxfc = -(66*u_mh - 85*uc_i + 23*uc_ip1 - 4*uc_ip2)/(18*dx);
    fh_imh =u_mh^2/2 - nu*uxfc ;%cell centered
    pretty(simplify(expand(fh_imh)))

    R_i1 = (fh_ciph-fh_imh)/dx;
    pretty(simplify(expand( R_i1 )))



%% Scheme consistent Neumann

    % Face centered
    uf0_3 = uf_ip1 - ux*dx;
    pretty(simplify(expand( uf0_3-uf_i )))

    % Boundary Flux
    fh_imh = (uf_ip1^2 + uf0_3^2)/4 - nu*(uf_ip1-uf0_3)/dx;
%     pretty(simplify(expand( fh_imh )))
    
    
    % Cell Centered
    uc0_3 = uc_i - ux_mh*dx;
%     pretty(simplify(expand( uc0_3-uc_im1 )))

    % Boundary Flux
    fh_imh = (uc_i^2 + uc0_3^2)/4 - nu*(uc_i-uc0_3)/dx;
%     pretty(simplify(expand( fh_imh )))

    % Residual for first cell
    R_i1 = (fh_ciph-fh_imh)/dx;
%     pretty(simplify(expand( R_i1 )))
    
    

    uf0_4 = (-12*ux*dx + 9*uf_ip1 + 3*uf_ip2 - uf_ip3)/(11);
%     pretty(expand( uf0_4-uf_i ))


    % Boundary Flux
    fh_imh = (uf_ip1^2 + uf0_4^2)/4 - nu*(uf_ip1-uf0_4)/dx;
%     pretty(simplify(expand( fh_imh )))
    
% Cell Centered


    uf0_4 = (-12*ux_mh*dx + 9*uc_i + 3*uc_ip1 - uc_ip2)/(11);
%     pretty(simplify(expand( uf0_4-uc_im1 )))

    ux_err = (uc_i-uf0_4)/dx - ux_mh;
%     pretty(simplify(expand(ux_err)));

    % Boundary Flux
    fh_imh = (uc_i^2 + uf0_4^2)/4 - nu*(uc_i-uf0_4)/dx;
%     fh_imh = (uc_i^2 + uf0_4^2)/4 - nu*ux_mh;
%     pretty(simplify(expand( fh_imh )))
    
    R_i1 = (fh_ciph-fh_imh)/dx;
%     pretty(expand( R_i1 ))

%% Scheme inconsistent Neumann

    % Face centered
    uf0_3 = -( 2*ux*dx - 7*uf_ip1 + uf_ip2 ) / 6;
    pretty(simplify(expand( uf0_3 )))

    % Boundary Flux
    fh_imh = uf0_3^2/2 - nu*ux;
    pretty(simplify(expand( fh_imh )))
    
    
    % Cell Centered
    uc0_3 = -( 2*ux_mh*dx - 7*uc_i + uc_ip1 ) / 6;
    pretty(simplify(expand( uc0_3-u_mh )))

    % Boundary Flux
    fh_imh = uc0_3^2/2 - nu*ux_mh;
    pretty(simplify(expand( fh_imh )))

    % Residual for first cell
    R_i1 = (fh_ciph-fh_imh)/dx;
    pretty(simplify(expand( R_i1 )))
    
    
    %======================================================================
    % 4th order extrapolation
    uf0_4 = -(18*ux*dx - 85*uf_ip1 + 23*uf_ip2 - 4*uf_ip3)/(66);
    pretty(expand( uf0_4 ))


    % Boundary Flux
    fh_imh = uf0_4^2/2 - nu*ux;
    pretty(simplify(expand( fh_imh )))
    
% Cell Centered

    uf0_4 = -(18*ux_mh*dx - 85*uc_i + 23*uc_ip1 - 4*uc_ip2)/(66);
    pretty(simplify(expand( uf0_4-u_mh )))

    % Boundary Flux
    fh_imh = uf0_4^2/2 - nu*ux_mh;
    pretty(simplify(expand( fh_imh )))
    
    R_i1 = (fh_ciph-fh_imh)/dx;
    pretty(expand( R_i1 ))



%% 

f1 = u + dx/2*ux + dx^2/8*uxx + dx^3/48*uxxx + dx^4*uxxxx/384 + dx^5*uxxxxx/3840;
f2 = u - dx/2*ux + dx^2/8*uxx - dx^3/48*uxxx + dx^4*uxxxx/384 - dx^5*uxxxxx/3840;
pretty( simplify( expand((f1^2+f2^2)/4 )) )
teface = simplify( expand( (f1^2+f2^2)/4) );

% teface = u^2/2 + dx^2/8*( ux^2 + u*uxx) + dx^4*(uxx^2/128 + ux*uxxx/96 + u*uxxxx/384);
% pretty(teface);

fc      = f     + dx/2*fx    + dx^2/8*fxx   + dx^3/48*fxxx  + dx^4*fxxxx/384;
fxc     = fx    + dx/2*fxx   + dx^2/8*fxxx  + dx^3/48*fxxxx + dx^4*fxxxxx/384;
fxxc    = fxx   + dx/2*fxxx  + dx^2/8*fxxxx + dx^3/48*fxxxxx;
fxxxc   = fxxx  + dx/2*fxxxx + dx^2/8*fxxxxx;
fxxxxc  = fxxxx + dx/2*fxxxxx;
fxxxxxc = fxxxxx;


tecellcenterp = subs( teface, {u, ux, uxx, uxxx, uxxxx, uxxxx, uxxxxx},{fc, fxc, fxxc, fxxxc, fxxxxc, fxxxxc, fxxxxxc});


pretty( simplify( expand( tecellcenterp ) ) )


% tecellcentern = subs( tecellcenterp, dx, -dx);
fc      = f     - dx/2*fx     + dx^2/8*fxx   - dx^3/48*fxxx  + dx^4*fxxxx/384;
fxc     = fx    - dx/2*fxx    + dx^2/8*fxxx  - dx^3/48*fxxxx + dx^4*fxxxxx/384;
fxxc    = fxx   - dx/2*fxxx   + dx^2/8*fxxxx - dx^3/48*fxxxxx;
fxxxc   = fxxx  - dx/2*fxxxx  + dx^2/8*fxxxxx;
fxxxxc  = fxxxx - dx/2*fxxxxx;
fxxxxxc = fxxxxx;


tecellcentern = subs( teface, {u, ux, uxx, uxxx, uxxxx, uxxxx, uxxxxx}, {fc, fxc, fxxc, fxxxc, fxxxxc, fxxxxc, fxxxxxc});



pretty( simplify( expand( tecellcentern ) ) )

TEinterior = (tecellcenterp-tecellcentern)/dx;

pretty(simplify(expand (TEinterior) ))


f1 = u + dx*ux + dx^2/2*uxx + dx^3/6*uxxx + dx^4*uxxxx/24 + dx^5*uxxxxx/120;
f2 = u - dx*ux + dx^2/2*uxx - dx^3/6*uxxx + dx^4*uxxxx/24 - dx^5*uxxxxx/120;

TEinterior_direct = (f1^2-f2^2)/(4*dx);

% f1 = f + dx*fx + dx^2/2*fxx + dx^3/6*fxxx + dx^4*fxxxx/24;
% f2 = f - dx*fx + dx^2/2*fxx - dx^3/6*fxxx + dx^4*fxxxx/24;
% TEdirectp = (f1^2+f^2)/4;
% TEdirectn = (f^2+f2^2)/4;

% pretty(expand(TEdirectp))
% pretty(expand(TEdirectn))
% pretty(expand( (TEdirectp-TEdirectn)/(2*dx) ))
pretty(expand(TEinterior_direct))
% 



% FV terms (face centered)
% Fh = [(u*ux)_(i+1/2) + (u*ux)_(i-1/2)]/4
fiph = u + ux/2*dx + uxx/6*dx^2 + uxxx/24*dx^3+uxxxx/120*dx^4+1/720*uxxxxx*dx^5;
fimh = u - ux/2*dx + uxx/6*dx^2 - uxxx/24*dx^3+uxxxx/120*dx^4-1/720*uxxxxx*dx^5;
pretty(expand( (fiph^2+fimh^2)/4 ))
pretty(expand( (fiph-fimh)/(dx) ))


% FV terms (cell centered face flux)
fiph = u + ux*dx + uxx*13/24*dx^2 + uxxx*5/24*dx^3+uxxxx*121/1920*dx^4+uxxxxx*91/5760*dx^5;
fimh = u + uxx/24*dx^2 + uxxxx/1920*dx^4;
pretty(expand(fiph^2))
pretty(expand(fimh^2))
FhR=(fiph^2+fimh^2)/4;
pretty(expand( FhR ))
pretty(expand( (fiph-fimh)/dx ))

fiph = u - ux*dx + uxx*13/24*dx^2 - uxxx*5/24*dx^3+uxxxx*121/1920*dx^4-uxxxxx*91/5760*dx^5;
fimh = u + uxx/24*dx^2 + uxxxx/1920*dx^4;
pretty(expand(fiph^2))
pretty(expand(fimh^2))
FhL=(fiph^2+fimh^2)/4;
pretty(expand( FhL ))

FhC = (FhR-FhL)/(dx);
pretty(expand(FhC))

fiph = u + ux*dx + uxx*13/24*dx^2 + uxxx*5/24*dx^3+uxxxx*121/1920*dx^4+uxxxxx*91/5760*dx^5 + uxxxxxx*1093/322560*dx^6;
fi   = u + uxx/24*dx^2 + uxxxx/1920*dx^4 + uxxxxxx/322560*dx^6;
fimh = u - ux*dx + uxx*13/24*dx^2 - uxxx*5/24*dx^3+uxxxx*121/1920*dx^4-uxxxxx*91/5760*dx^5 + uxxxxxx*1093/322560*dx^6;
pretty(expand( - nu*(fiph-2*fi+fimh)/(dx^2) ))


fprintf('Full residual for interior scheme:\n\n')
pretty(expand( (fiph^2-fimh^2)/(4*dx) - nu*(fiph-2*fi+fimh)/(dx^2) ))

fprintf('Full residual for dirichlet BC scheme consistent formulation:\n\n')
utbar = ut + uxxt/24*dx^2 + uxxxxt/1920*dx^4 + uxxxxxxt/322560*dx^6
% pretty(expand( (fiph^2-utbar^2)/(4*dx) - nu*(fiph-2*fi+utbar)/(dx^2) ))
% pretty(expand( (ut^2+fi^2)/4 ))
% pretty(expand( -nu*(fi-utbar)/dx ))

%centered at the face
fface = u
fi = u + ux/2*dx + uxx/6*dx^2 + uxxx/24*dx^3+uxxxx/120*dx^4+1/720*uxxxxx*dx^5;
fiph = u + 3/2*ux*dx + 7/6*uxx*dx^2 + 5/8*uxxx*dx^3 + 31/120*uxxxx*dx^4 + 7/80*uxxxxx*dx^5;
pretty(expand( ( -( 6*fface - 7*fi + fiph)/(2*dx) ) ) )

% cell centered
fface = u - dx/2*ux + dx^2/8*uxx -dx^3/48*uxxx + dx^4/384*uxxxx - dx^5/3840*uxxxxx + dx^6/46080*uxxxxxx
fi   = u + uxx/24*dx^2 + uxxxx/1920*dx^4 + uxxxxxx/322560*dx^6;
fiph = u + ux*dx + uxx*13/24*dx^2 + uxxx*5/24*dx^3+uxxxx*121/1920*dx^4+uxxxxx*91/5760*dx^5 + uxxxxxx*1093/322560*dx^6;
pretty(expand( (-( 6*fface - 7*fi + fiph)/(2*dx) ) ) )

% fiph = u + ux/2*dx + uxx/6*dx^2 + uxxx/24*dx^3+uxxxx/120*dx^4+1/720*uxxxxx*dx^5;
% fip3h = u + 3/2*ux*dx + 7/6*uxx*dx^2 + 5/8*uxxx*dx^3 + 31/120*uxxxx*dx^4 + 7/80*uxxxxx*dx^5;
% pretty(expand( -( 6*u - 7*fiph + fip3h)/(2*dx) ) )

fprintf('Dirichlet scheme inconsistent\n')
pretty(expand(fface^2/2))
pretty(expand( fface^2/2 - nu*(-( 6*fface - 7*fi + fiph)/(2*dx)) ) )

Fl = fface^2/2 - nu*(-( 6*fface - 7*fi + fiph)/(2*dx))
Fr = (fiph^2+fi^2)/4 - nu*(fiph-fi)/(dx);
pretty(expand(Fr))

pretty(expand( (Fr-Fl)/dx ))

% 
% fprintf('Neumann scheme inconsistent\n  Centered at the face\n')
% fiph = u + ux/2*dx + uxx/6*dx^2 + uxxx/24*dx^3+uxxxx/120*dx^4+1/720*uxxxxx*dx^5;
% fip3h = u + 3/2*ux*dx + 7/6*uxx*dx^2 + 5/8*uxxx*dx^3 + 31/120*uxxxx*dx^4 + 7/80*uxxxxx*dx^5;
% pretty(expand( -1/3*dx*ux + 7/6*fiph - 1/6*fip3h ) )
% 
% fprintf('  Centered at the cell center\n')
% dfface =  ux - dx/2*uxx + dx^2/8*uxxx -dx^3/48*uxxxx + dx^4/384*uxxxxx - dx^5/3840*uxxxxxx;
% fi   = u + uxx/24*dx^2 + uxxxx/1920*dx^4 + uxxxxxx/322560*dx^6;
% fiph = u + ux*dx + uxx*13/24*dx^2 + uxxx*5/24*dx^3+uxxxx*121/1920*dx^4+uxxxxx*91/5760*dx^5 + uxxxxxx*1093/322560*dx^6;
% fface = (-2*dx*dfface + 7*fi - fiph)/6;
% pretty(expand(fface));
% pretty(expand( fface^2/2 - nu*dfface  ))