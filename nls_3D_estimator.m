function [np_est, residual_norm2] = nls_3D_estimator(r,a,initial_estimate, Niter, w, x1,x2,weights,np)

np_est = initial_estimate;
% np_est = np+randn(3,1);
Na = size(a,2);
for iter = 1:Niter
%     np_est(1) = np(1);
%     X1(1) = np_est(1) + wh/2;
%     X2(1) = X1(1);
    H = zeros(Na,3);
    p = zeros(Na,1);
    xn = np_est(1);
    yn = np_est(2);
    zn = np_est(3);
    zb = zn-w/2;
    for aidx = 1:Na
            [~,pl,dp_dxn,dp_dyn,dp_dzn] = get_diffraction_coord_fermat(a(1,aidx),a(2,aidx),a(3,aidx),x1,x2,xn,yn,zn,zb,0.001);
            H(aidx,1) = dp_dxn;
            H(aidx,2) = dp_dyn;
            H(aidx,3) = dp_dzn;
            p(aidx) = pl;
    end
    weights = eye(Na);
    Htw = transpose(H)*weights;
    HtwH = Htw*H;
    
    % H_pseudo = (transpose(H)*weights*H)^-1*transpose(H)*weights;
    H_pseudo = (HtwH)^-1*transpose(H)*weights;
    residual_2 = (r-p);
    residual_norm2 = sqrt(sum(residual_2.^2));
    if cond(HtwH)>1000 && det(HtwH) < 0.1
        np_est = np_est + 1/(trace(HtwH))*Htw*(residual_2);
        
    else
    np_est = np_est + H_pseudo*(residual_2);
    end
    
    % iter
    % cond(H_pseudo)]
    % [initial_estimate,np_est,np]
    
end

if sqrt(sum((np-np_est).^2)) > 1e3
    disp('error')
end
end