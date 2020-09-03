function Pcssm=CSSM(X,A,grid,K)
c = 3e8;
L = 10;             % Number of sensors
u = [12 40 48];
nS = length(u);     % Number of sources
fs = 1; Ts = 1/fs;
fc = 0.3; B = 0.2*ones(1,nS);
d = 0.5*c/(max(fc)+max(B)/2); lam = c/(max(fc)+max(B)/2);
d1 = 3*d; d2 = 4*d; 
ang = (0:4)*2*pi/5;
elev = pi/2;
AG = zeros(3,L);
AG(:,1:2:L) = d1*[cos(ang')*sin(elev) sin(ang')*sin(elev) cos(elev)*ones(size(ang'))]';
AG(:,2:2:L) = d2*[cos(ang')*sin(elev) sin(ang')*sin(elev) cos(elev)*ones(size(ang'))]';
beamwidth = lam/(2*d2)*180/pi;

Restta = .125/6;
Tta = 1:Restta:60; 

F = 256;            % FFT points

fL = fc - max(B)/2; fU = fc + max(B)/2;
nfc = round(fc*F); fc = nfc/F;
nfL = round(fL*F); nfU = round(fU*F);
fL = nfL/F; fU = nfU/F;
f = (fL:fs/F:fU);
lf = length(f);

B_use = 0.2;
nbins = round(B_use*F);
if mod(nbins,2) == 0
    nbins = nbins + 1;
end

lr = nfc - floor(nbins/2); ur = nfc + floor(nbins/2);

J = 10;             % Number of frequency-domain snapshots
t = Ts:Ts:J;  lt = length(t);

        % Signal generation
        at = sqrt(randn(nS,lf).^2+randn(nS,lf).^2);
        for kkk = 1:lf
            for k = 1:nS
                s(k,:,kkk) = at(k,kkk).*exp(j*( 2*pi*( f(kkk))*t +(2*rand(1,length(t))-1)*pi));
                a(:,k,kkk) = spv(AG,[0 u(k)], f(kkk), c);
            end

            for ll = 1:L
                n(ll,:,kkk) = randn(1,lt)+j*randn(1,lt);
                n(ll,:,kkk) = n(ll,:,kkk)/sqrt(var(n(ll,:,kkk))); %% Mistake cleared - replaced norm with var
            end
        end
     
        
        for kkk = 1:lf
            x(:,:,kkk) = a(:,:,kkk)*s(:,:,kkk)*sqrt(snr)+n(:,:,kkk);
        end
        Xf = randn(J,L,F);
        for k = 1:L
            Xf(:,k,nfL:nfU) = x(k,:,:);
        end

        Xff = [];
        for k = lr:ur
            Xff = [Xff Xf(:,:,k).'];
        end
        R_f = (Xff-mean(Xff,2)*ones(1,nbins*J))*(Xff-mean(Xff,2)*ones(1,nbins*J))'/(nbins*J);
        [Uaa Vaa] = eigs(R_f,L);
        
        % Initial estimates for CSSM
        for tta = 1:length(Tta)
            Pmvdr(tta) = 10*log10(1./abs(diag(sv(:,tta)'*inv(R_f)*sv(:,tta))).^1);
        end
        
        bb = findpeaks1(Pmvdr);
        azi = [];
        for kb = 1:length(bb)
            azi = [azi Tta(bb(kb)) Tta(bb(kb))-0.125*beamwidth Tta(bb(kb))+0.125*beamwidth];
        end
        zia = [cosd(azi);sind(azi);zeros(1,length(azi))];
        
        % CRLB evaluation and CSSM, ISSM evaluation
        temp = zeros(L); Rn  = zeros(L); Ry = zeros(L);
        for kkk = 1:lf
            Rx = (x(:,:,kkk)-mean(x(:,:,kkk),2)*ones(1,J))*(x(:,:,kkk)-mean(x(:,:,kkk),2)*ones(1,J))'/J;
            iRx = inv(Rx);   
            % For cssm
            ksv = exp(-i*2*pi*f(kkk)/c*AG'*zia);
            ksv0 = exp(-i*2*pi*fc/c*AG'*zia);
            prod = ksv*ksv0';
            [Uf D Vf] = svd(prod);
            T = Vf*Uf';
            Ryf{kkk} = T*Rx*T';
            [aa bb] = eigs(Rx,L);
            db = diag(bb);
            sigv(kkk) = mean(db(nS+2:end));
            Ry = Ry + Ryf{kkk}/lf;
            Rn = Rn + sigv(kkk)*T*T';   
        end    
        % For cssm
        [Ucssm Vcssm] = eigs(inv(Rn)*Ry,L);
        Uncssm = Ucssm(:,nS+1:end);
        for tta = 1:length(Tta);
            P_cssm(tta) = 10*log10(1./abs(sv(:,tta)'*Uncssm*Uncssm'*sv(:,tta)));
        end