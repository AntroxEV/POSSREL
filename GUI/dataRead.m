function [Mqt] = dataRead(filename,filetype,zLayer)
% reading database and elaboration (uniform discretization and correction
% pore pressure)
load(filename)
anet=0.6;
indX=1;
dz=zLayer(2)-zLayer(1);
z1=zLayer(1);
z2=zLayer(end);
zz=0:dz:z2;
Mqt=[];
switch filetype
    case 0
        
        for hh=1:nID
            vect=unique(HOLEID);
            if isstring(vect)
                indx=categorical(HOLEID)==vect{hh};
            else
                indx=HOLEID==vect(hh);
            end
            MM=sortrows([STCNDPTH(indx) STCNRES(indx) STCNPWP1(indx)]);
            qc=MM(:,2).*1E6; %convert to Pa CONE RESISTANCE
            z=MM(:,1);
            pw=MM(:,3).*1E3; %convert to Pa PORE PRESSURE
            pw(end)=0; %correction
            if length(z) ~= length(unique(z))
                % z=z(1):z(2)-z(1):z(end);
                [au,ia] = unique(z,'stable');
                z=z(ia);
                qc=qc(ia);
                pw=pw(ia);
            end
            
            if z(end)>z2 %discarge short CPT
                ind=zz<=z1;
                if ~isempty(ind)
                    qt=qc+pw*(1-anet); %correction (corrected cone resistance)
                    % create homogeneous optimization
                    qti=interp1(z,qt,zz,"nearest");
                    qti(isnan(qti))=0;
                    qti(ind)=0;
                    qti(qti<0)=0;
                    Mqt(:,indX)=qti; %saving
                    % plot(qti,-zd,'Color',[.7 .7 .7]);
                    indX=indX+1;
                end
            end
        end

    case 1
        NIDcases=unique(MDATASETS(2:end,1));
        nID=length(NIDcases);
        for hh=1:nID
            indx=MDATASETS(:,1)==NIDcases(hh);
            qti=MDATASETS(indx,3).*1E6; %convert to Pa CONE RESISTANCE
            z=MDATASETS(indx,2);
            if length(z) ~= length(unique(z))
                % z=z(1):z(2)-z(1):z(end);
                [au,ia] = unique(z,'stable');
                z=z(ia);
                qti=qti(ia);
            end

            if z(end)>z2 %discarge short CPT
                ind=zz<=z1;
                if ~isempty(ind)
                    qti=interp1(z,qti,zz,"nearest");
                    qti(isnan(qti))=0;
                    qti(ind)=0;
                    qti(qti<0)=0;
                    Mqt(:,indX)=qti; %saving
                    indX=indX+1;
                end
            end
        end
end
if isempty(Mqt)
    Mqt=[];
end