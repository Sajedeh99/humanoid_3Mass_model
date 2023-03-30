function [qpresult, Opt_Vector]=controller2(t, T, Lnom, Wnom, L_min, L_max, W_min, W_max, T_min, T_max,...
    bx_nom, by_nom, w0, zeta_mea_x, zeta_mea_y, uTx_ref, uTy_ref, zeta_err_x, zeta_err_y, PcZMP_y, PcZMP_x)

Wux = 1;
Wuy = 1;
WT = 0.05;
Wbx = 5;
Wby = 5;

H=[ 2*Wux    0        0        0        0
    0        2*Wbx    0        0        0
    0        0        2*Wuy    0        0
    0        0        0        2*Wby    0
    0        0        0        0        2*WT];
f=[];

A=[
     1     0     0     0     0
    -1     0     0     0     0
     0     0     1     0     0
     0     0    -1     0     0
     0     0     0     0     1
     0     0     0     0    -1
    ];

B=[
     L_max - uTx_ref
   -(L_min - uTx_ref)
     W_max - uTy_ref
   -(W_min - uTy_ref)
     exp(w0*T_max) - exp(w0*T)
   -(exp(w0*T_min) - exp(w0*T))
    ];

Aeq=[
    1 1 0 0 -(zeta_mea_x(2) - PcZMP_x)*exp(-w0*t);
    0 0 1 1 -(zeta_mea_y(2) - PcZMP_y)*exp(-w0*t)
    ];
Beq=[(zeta_err_x(2)-PcZMP_x)*exp(w0*(T-t))+PcZMP_x,(zeta_err_y(2)-PcZMP_y)*exp(w0*(T-t))+PcZMP_y]';

qpresult = quadprog(H, f, A, B, Aeq, Beq);
Opt_Vector(1,1) = qpresult(1) + uTx_ref;
Opt_Vector(2,1) = qpresult(3) + uTy_ref;
Opt_Vector(3,1) = qpresult(5) + exp(w0*T);
Opt_Vector(4,1) = qpresult(2) + bx_nom;
Opt_Vector(5,1) = qpresult(4) + by_nom;
end