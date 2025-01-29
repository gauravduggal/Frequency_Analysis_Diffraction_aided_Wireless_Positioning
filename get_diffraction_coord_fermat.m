function [Qe,p,dp_dxn,dp_dyn,dp_dzn] = get_diffraction_coord_fermat(xa,ya,za,x1,x2,xn,yn,zn,zb,delta)
%Secret Sauce for Diffraction path model: Validated using numeric Fermat's
% principle of least time. % Refer to https://arxiv.org/pdf/2409.02832 for
% derivation

a = (x1-x2)^2*((yn^2-ya^2)+(zn-za)*(zn+za-2*zb));
b = 2*(x1-x2)*((x2-xa)*((zb-zn)^2+yn^2)-(x2-xn)*((zb-za)^2+ya^2));
c = (x2-xa)^2*((zb-zn)^2+yn^2)-(x2-xn)^2*((zb-za)^2+ya^2);
D = sqrt(b^2-4*a*c);
qy = 0;
qz = zb;

root_flag = -1;
degenerate_flag = 0;
lambda = (-b-D)/(2*a);
if (b^2-4*a*c)<delta
    lambda = -b/(2*a);
    degenerate_flag = 1;
end

if abs(a)<delta
    lambda = -c/b;
    degenerate_flag = 2;
    if abs(b)< delta
        lambda = 0;
        degenerate_flag = 3;
    end
end

%if two roots
if degenerate_flag == 0
    lambda_pos_root = (-b+D)/(2*a);
    qx_pos_root = lambda_pos_root*x1+(1-lambda_pos_root)*x2;
    Qe_pos_root = [qx_pos_root;qy;qz];

    lambda_neg_root = (-b-D)/(2*a);
    qx_neg_root = lambda_neg_root*x1+(1-lambda_neg_root)*x2;
    Qe_neg_root = [qx_neg_root;qy;qz];

    [p_pos_root, IPL_pos_root, OPL_pos_root] = get_diff_length(Qe_pos_root,[xa;ya;za],[xn;yn;zn]);
    [p_neg_root, IPL_neg_root, OPL_neg_root] = get_diff_length(Qe_neg_root,[xa;ya;za],[xn;yn;zn]);

    if p_pos_root < p_neg_root
        root_flag = 1;
        Qe = Qe_pos_root;
        IPL = IPL_pos_root;
        OPL = OPL_pos_root;
        p = p_pos_root;
        qx = qx_pos_root;
    else
        root_flag = -1;
        Qe = Qe_neg_root;
        IPL = IPL_neg_root;
        OPL = OPL_neg_root;
        p = p_neg_root;
        qx = qx_neg_root;
    end
else
    qx = lambda*x1+(1-lambda)*x2;
    Qe = [qx;qy;qz];
    [p,IPL,OPL] = get_diff_length(Qe,[xa;ya;za],[xn;yn;zn]);
end



dc_dzn = -2*(x2-xn)^2*(zb-za);
dc_dyn = 2*yn*(x2-xa)^2;
dc_dxn = 2*(x2-xn)*((zb-za)^2+ya^2);

db_dzn = -4*(x1-x2)*(x2-xn)*(zb-za);
db_dyn = 4*yn*(x1-x2)*(x2-xa);
db_dxn = 2*(x1-x2)*((zb-za)^2+ya^2);

da_dxn = 0;
da_dyn = 2*yn*(x1-x2)^2;
da_dzn = (x1-x2)^2*(2*(za-zb));

if root_flag == 0
    disp('error')
end

if degenerate_flag == 0
    dqx_dzn = (x1-x2)/(2*a)*(da_dzn*(b-root_flag*D)/a+(-db_dzn+root_flag*(b*db_dzn-2*c*da_dzn-2*a*dc_dzn)/D));
    dqx_dyn = (x1-x2)/(2*a)*(da_dyn*(b-root_flag*D)/a+(-db_dyn+root_flag*(b*db_dyn-2*c*da_dyn-2*a*dc_dyn)/D));
    dqx_dxn = (x1-x2)/(2*a)*(da_dxn*(b-root_flag*D)/a+(-db_dxn+root_flag*(b*db_dxn-2*c*da_dxn-2*a*dc_dxn)/D));
elseif degenerate_flag == 1
    dqx_dzn = (x1-x2)*(-db_dzn*a+b*da_dzn)/(2*a);
    dqx_dyn = (x1-x2)*(-db_dyn*a+b*da_dyn)/(2*a);
    dqx_dxn = (x1-x2)*(-db_dxn*a+b*da_dxn)/(2*a);
elseif degenerate_flag == 2
    dqx_dzn = (x1-x2)*(-dc_dzn*b+c*db_dzn)/(b^2);
    dqx_dyn = (x1-x2)*(-dc_dyn*b+c*db_dyn)/(b^2);
    dqx_dxn = (x1-x2)*(-dc_dxn*b+c*db_dxn)/(b^2);
else
    dqx_dzn = 0;
    dqx_dyn = 0;
    dqx_dxn = 0;
end
dp_dzn = ((qx-xa)*dqx_dzn+(zb-za))/OPL - (xn-qx)*dqx_dzn/IPL;
dp_dyn = (qx-xa)*dqx_dyn/OPL + (-(xn-qx)*dqx_dyn+yn)/IPL;
dp_dxn = (qx-xa)*dqx_dxn/OPL + (xn-qx)*(1-dqx_dxn)/IPL;



end

function [p,IPL,OPL] = get_diff_length(Qe,Xa,Xn)
si = Qe-Xa;
sd = Xn-Qe;
OPL = sqrt(sum(si.^2));
IPL = sqrt(sum(sd.^2));
p = OPL + IPL;
end

