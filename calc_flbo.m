function [W,A] = calc_flbo(vertices, faces, options)

n = size(vertices,1);
m = size(faces,1);

% Reorient mesh faces if they are inconsistent

adjacency_matrix = sparse([faces(:,1); faces(:,2); faces(:,3)], ...
                         [faces(:,2); faces(:,3); faces(:,1)], ...
    	                 ones(3 * m, 1), ...
                         n, n, 3 * m);
if any(any(adjacency_matrix > 1))
    options.method = 'slow';
    warning('Inconsistent face orientation. The mesh will be reoriented.')
    faces = transpose(perform_faces_reorientation(vertices,faces,options));
end
clear adjacency_matrix

[U1, U2, D, N] = avg_diffusion_tensor_flbo(vertices, faces, options.alpha, options.curv_smooth, options.angle);

% Construct the Finsler-based stiffness matrix

W = sparse(n,n);

angles = zeros(size(faces));
for i=1:3
    i1 = mod(i-1,3)+1;
    i2 = mod(i  ,3)+1;
    i3 = mod(i+1,3)+1;
    pp = vertices(faces(:,i2),:) - vertices(faces(:,i1),:);
    qq = vertices(faces(:,i3),:) - vertices(faces(:,i1),:);
    pp = pp ./ repmat( max(sqrt(sum(pp.^2,2)),eps), [1 3] );
    qq = qq ./ repmat( max(sqrt(sum(qq.^2,2)),eps), [1 3] );
    angles(:,i1) = acos(sum(pp.*qq,2));
end



for i=1:3
    i1 = mod(i-1,3)+1;
    i2 = mod(i  ,3)+1;
    i3 = mod(i+1,3)+1;
    % we only focus on (i1,i1) and on (i1,i2)
    e1 = vertices(faces(:,i3),:)-vertices(faces(:,i2),:);
    e2 = vertices(faces(:,i1),:)-vertices(faces(:,i3),:);
    e1 = e1./repmat(sqrt(sum(e1.^2,2)),[1 3]);
    e2 = e2./repmat(sqrt(sum(e2.^2,2)),[1 3]);

    tau = options.tau; % Randers coefficient
    U1rep = repmat(U1,1,3);
    U2rep = repmat(U2,1,3);

    Nrep = repmat(N,1,3);
    Nelem = repelem(N,1,3);

    U1elem = repelem(U1,1,3);
    U2elem = repelem(U2,1,3);
    w = tau * U2; % Example of w for the (M,w) couple that defines Randers metric
    D(:,3) = 1;
    w_RiemmanProduct =    D(:,1).*(sum(w.*U1, 2).^2) +  D(:,2).*(sum(w .* U2, 2).^2) + D(:,3).*(sum(w.* N, 2).^2); % w_RiemmanProduct = < w, M.inverse() * w >
    eta = 1 - w_RiemmanProduct;

    UDUT = D(:,1) .* (U1elem .* U1rep) + D(:,2) .* (U2elem .* U2rep) + D(:,3) .* (Nelem .* Nrep); % Derive the shear matrix


    UDUTw(:,1) = sum(UDUT(:,1:3) .* w(:,1:3), 2);
    UDUTw(:,2) = sum(UDUT(:,4:6) .* w(:,1:3), 2);
    UDUTw(:,3) = sum(UDUT(:,7:9) .* w(:,1:3), 2);

    UDUTwelem = repelem(UDUTw, 1, 3);
    UDUTwrep = repmat(UDUTw, 1,3);

    Mstar = (UDUTwrep .* UDUTwelem + eta .* UDUT) ./ (eta .* eta); % Derive M*, dual of M
    wstar = - UDUTw ./ eta; % Derive w*, dual of w

    wstarrep = repmat(wstar,1,3);
    wstarelem = repelem(wstar,1,3);


    wProdw = wstarrep .* wstarelem;
    DFinsler = Mstar - wProdw; 
    e1elem = repelem(e1,1,3);
    e1rep = repmat(e1,1,3);
    e2elem = repelem(e2,1,3);
    e2rep = repmat(e2,1,3);

    DFinsler(abs(DFinsler)<1e-10) = 0;
    DFinsler(abs(DFinsler-1)<1e-10) = 1;


    factore = -(1/2) * ...
        sum(e1elem .* (DFinsler .* e2rep), 2) ./...
        sin(angles(:,i3));

    factord = -(1/2) * ...
        sum(e1elem .* (DFinsler .* e1rep), 2) .*...
        (cot(angles(:,i2))+cot(angles(:,i3)));

    W = W + sparse([faces(:,i1);  faces(:,i2); faces(:,i1)],...
                   [faces(:,i2);  faces(:,i1); faces(:,i1)],...
                   [factore;      factore;     factord    ],...
                        n, n);

end
% Construct the mass matrix
S_tri = zeros(m,1);
for k=1:m
    e1 = vertices(faces(k,3),:) - vertices(faces(k,1),:);
    e2 = vertices(faces(k,2),:) - vertices(faces(k,1),:);
    S_tri(k) = 0.5*norm(cross(e1,e2));
end
A = zeros(n,1);
for i=1:m
    A(faces(i,1)) = A(faces(i,1)) + S_tri(i)/3;
    A(faces(i,2)) = A(faces(i,2)) + S_tri(i)/3;
    A(faces(i,3)) = A(faces(i,3)) + S_tri(i)/3;
end
A = sparse(1:n,1:n,A);

end
