%% Par�metros de Sistema

clear;
NumSub=52;  %Quantidade de subportadoras do sistema OFDM
Qsimb=100;   %Quantidade de s�mbolos � ser transmitido
Nfft = 2^(nextpow2(NumSub));%Total de subportadoras usadas e n�o-usadas
Pcp=25; %Quantidade de prefixo ciclico ou intervalo de guarda
Nbg=Nfft-NumSub; % Banda Virtual (ZeroPading)



%% Configura��o
TipoMod=1; %Tipo de modula��o %1->16-QAM  % 2->BPSK
TipoCanal=2; %Tipo de canal %1->AWGN 2->Rayleigh
IntGuard=1; %0-> N�o tiver Int. Guard %1 -> Se tiver Int. Guarda
            Ogi=0; %1 -> Intervalo de Guarda (Zero Padding)
            Opc=1; %1 -> Prefixo C�clico

            
            
            
%% Variaveis de Erros
EbNoin = -10;
EbNofim = 35;
EbNo_dB = EbNoin:5:EbNofim; % Varia��o de Eb/N0 dB
    
Ncp = (Pcp/100)*Nfft; %Tamanho Prefixo C�clico


Qrep=10; % N�mero de Repeti��es

%Erros

snr=zeros(1,length(EbNo_dB));
Erros=zeros(1,length(EbNo_dB));
ErrosG=zeros(1,length(EbNo_dB));


for i=1:length(EbNo_dB)
    for k=1:Qrep
        %Iniciar as Vari�veis
        vtBg=(1:NumSub/2);
        vtBgf=(floor(NumSub/2)+1:NumSub);
        vtFft=(1:Nfft);
        vtRBg=(2:(floor(NumSub/2)+1));
        vtRBgf=(floor(NumSub/2)+(Nbg+1):Nfft);
        vtCp=(Nfft-(Ncp-1):Nfft);
        vtRCp=(Ncp+1:(NumSub+Nbg+Ncp));

        xfreq=[];
        xtime=[];
        yfreq=[];
        yfreqG=[];
        ytime=[];
        ytimeG=[];
        Tx=[];
        TxG=[];
        Rx=[]; 
        RxG=[];
    
        
%% Gerador de Bits e Modula��o

        if TipoMod==1 % Modula��o 16-QAM
            
            Mod = 4; %Tipo de modula��o (4 bits/simbolos)
            Qbits= Mod*NumSub*Qsimb; %Quantidade de bits
            a=randi([0 1],1,Qbits); %Gerando Bits aleat�rios
            ainf=bi2de(transpose(reshape(a,Mod,length(a)/Mod)),'left-msb');
        
            mod=qammod(transpose(ainf),16,0,'gray')/sqrt(10); %Modula��o 16-Qam

        elseif TipoMod==2 % Modula��o BPSK
            Mod = 1; %Tipo de modula��o (4 bits/simbolos)
            Qbits= Mod*NumSub*Qsimb;
            mod=2*round(rand(1,Qbits))-1; %Gerando Bits em modula��o BPSK
        end
        
%% Banda de Guarda (Banda Virtual)
        for m=1:Qsimb

             xfreq=[xfreq [0 mod(vtBg) zeros(1,Nbg-1) mod(vtBgf)]];
             vtBg=vtBg+NumSub;
             vtBgf=vtBgf+NumSub;
        end
        
        
%% IFFT
        for m=1:Qsimb
            xtime = [ xtime ifft(xfreq(vtFft))];

            vtFft=vtFft+Nfft;
        end

        vtFft=(1:Nfft); % Reiniciar vari�vel
        
        %% Adicionar Prefixo C�clico/Intervalo de Guarda
        if IntGuard==1 && Opc==1 && Ogi==0  % Se for Prefixo C�clico
    
            % Prefixo C�clico
            for m=1:Qsimb
                TxG = [TxG [xtime(vtCp) xtime(vtFft)]];

                vtCp=vtCp+Nfft;
                vtFft=vtFft+Nfft;
            end
            
        elseif IntGuard==1 && Opc==0 && Ogi==1  %Se for Intervalo de Guarda
            %Preencher com Zeros
             for m=1:Qsimb
                TxG = [TxG [zeros(1,Ncp) xtime(vtFft)]];

                vtFft=vtFft+Nfft;
            end
            

        end
        Tx = xtime; % Sem Intervalo de Guarda e nem prefixo C�clico
        
  %% Modelagem do Canal   
        
 %% Modelagem do Canal   
        
       
       snr(i)=EbNo_dB(i) + 10*log10(Mod) + 10*log10(NumSub/Nfft);
       
        if IntGuard==1 % Se tiver Intervalo de Guarda(Zero Padding) ou Prefixo C�clico
        Eb=(TxG*TxG')/length(TxG);
        
        %Rayleigh
                    
        if TipoCanal==2 % Se tiver Espalhaento Rayleigh
      
                hG = (1/sqrt(2))*(randn(1,1))+1i*randn(1,1);
                TxG = conv(TxG,hG);
%                 TxG =  TxG.*rayG;     
                   
        end
            
            No=Eb*10.^(-snr(i)/10);
            Rx_noiseG = TxG+(sqrt(No/2))*((randn(1,length(TxG)) + 1j*randn(1,length(TxG))));     
        end
        
        %AWGN
        Eb=(Tx*Tx')/length(Tx);
        
        if TipoCanal==2 % Se tiver Espalhaento Rayleigh
        h = (1/sqrt(2))*(randn(1,1))+1i*randn(1,1);
        Tx = conv(Tx,h); 
        end
        
        No=Eb*10.^(-snr(i)/10);
        Rx_noise = Tx+(sqrt(No/2))*((randn(1,length(Tx)) + 1j*randn(1,length(Tx))));
        
         
        
         
         
%% Removendo Prefixo C�clico/Intervalo de Guarda
         Rx = Rx_noise;
         
         if IntGuard==1 
             for m=1:Qsimb
                RxG = [RxG Rx_noiseG(vtRCp)];

               vtRCp=vtRCp+(Nfft+Ncp);
  
             end
         end
         
         
%% FFT
         
         vtFft=(1:Nfft); %Reiniciando Vari�vel
         
        %FFT sem Banda de Guarda
         for m=1:Qsimb
            yfreq = [ yfreq fft(Rx(vtFft))];
            vtFft=vtFft+Nfft;
         end
        
        
        if IntGuard==1
            %FFT com Banda de Guarda
            vtFft=(1:Nfft); %Reiniciando Vari�vel
            
            for m=1:Qsimb
                yfreqG = [ yfreqG fft(RxG(vtFft))];
                vtFft=vtFft+Nfft;
            end
        end
        
        
        
%% Retirada da Banda de Guarda (Banda Virtual)
         vtRBg=(2:(floor(NumSub/2)+1)); %Reiniciando Vari�vel
        vtRBgf=(floor(NumSub/2)+(Nbg+1):Nfft);
        
        for m=1:Qsimb
             ytime=[ytime [yfreq(vtRBg) yfreq(vtRBgf)]];
             vtRBg=vtRBg+Nfft;
             vtRBgf=vtRBgf+Nfft;
        end
        
        if IntGuard==1
             vtRBg=(2:(floor(NumSub/2)+1)); %Reiniciando Vari�vel
            vtRBgf=(floor(NumSub/2)+(Nbg+1):Nfft);
            for m=1:Qsimb
                 ytimeG=[ytimeG [yfreqG(vtRBg) yfreqG(vtRBgf)]];
                 vtRBg=vtRBg+Nfft;
                 vtRBgf=vtRBgf+Nfft;
            end
        end
  
        %%Equaliza��o caso (Rayleigh)
         if TipoCanal==2 % Se tiver Espalhaento Rayleigh

            
            if IntGuard==1 %Se tiver intervalo de guarda
                
                vtRBg=(2:(floor(NumSub/2)+1)); %Reiniciando Vari�vel
                vtRBgf=(floor(NumSub/2)+(Nbg+1):Nfft);

                hFG= fft(hG,Nfft);
                hFG=[hFG(vtRBg) hFG(vtRBgf)];
                vtEq = (1:NumSub); 

                for m = 1:Qsimb                
                ytimeG(vtEq) = ytimeG(vtEq)./hFG; %equalizador
                vtEq = vtEq + NumSub; 
                end
            end
            
            vtRBg=(2:(floor(NumSub/2)+1)); %Reiniciando Vari�vel
            vtRBgf=(floor(NumSub/2)+(Nbg+1):Nfft);
            
            hF= fft(h,Nfft);
            hF=[hF(vtRBg) hF(vtRBgf)];
            vtEq = (1:NumSub); 
            
            for m = 1:Qsimb                
            ytime(vtEq) = ytime(vtEq)./hF; %equalizador
            vtEq = vtEq + NumSub;
            end           
         end
        
        
%% Demodula��o
        
        if TipoMod==1  %Demodula��o 16-QAM
            
          b=qamdemod(ytime*sqrt(10),16,0,'gray');
          binf2=de2bi(b,'left-msb');       
          binf=reshape(transpose(binf2),[1,length(a)]);
          
          if IntGuard==1 %Se tiver Int. Guarda
              bG=qamdemod(ytimeG*sqrt(10),16,0,'gray');
              binf2G=de2bi(bG,'left-msb');       
              binfG=reshape(transpose(binf2G),[1,length(a)]);
              
              %Erros
              errG=sum(a~=binfG);
              ErrosG(i)=ErrosG(i)+errG;
          end
          
            %Erros
            err=sum(a~=binf);
            Erros(i)=Erros(i)+err;

        elseif TipoMod==2    %Demodula��o BPSK
            
            infOut2=real(ytime);
            infOut2(infOut2>0) = +1;
            infOut2(infOut2<0) = -1;
            infOut=infOut2;
            
            if IntGuard==1 %Se tiver Int. Guarda
                
                infOut2G=real(ytimeG);
                infOut2G(infOut2G>0) = +1;
                infOut2G(infOut2G<0) = -1;
                infOutG=infOut2G;
                
                %Erros
                errG=sum(a~=infOutG);
                ErrosG(i)=ErrosG(i)+errG;
            end
            
            %Erros
            err=sum(a~=infOut);
            Erros(i)=Erros(i)+err;
        end      
        
        
        
    end
end


%% C�lculo BER

EbNoLin=10.^(EbNo_dB/10); %Eb/N0 em escala linear
gamma_c=EbNoLin*Mod;

BerTeoricoBpskAwgn=(1/2)*erfc(sqrt(10.^(EbNo_dB/10))); % BER BPSK AWGN
berTeoricoBpskRay=0.5*(1-sqrt(EbNoLin./(1+EbNoLin))); %BER BPSK Rayleigh

BerTeoricoQamAWGN = (1/Mod)*3/2*erfc(sqrt(Mod*0.1*(10.^(EbNo_dB/10))));%BER 16-QAM AWGN
berTeoricoQam16Ray = 3/8 * ( 1 - sqrt(2/5*gamma_c/Mod./(1+2/5*gamma_c/Mod)) ); %BER 16-QAM Rayleigh

BerSimulado = Erros/(Qsimb)/(Mod)/(Qrep)/(NumSub); %BER SIMULADO

if IntGuard==1 %Se tiver Int. Guarda
    BerSimuladoG = ErrosG/(Qsimb)/(Mod)/(Qrep)/(NumSub); %BER SIMULADO
end

%% Plot

figure(1);
semilogy(EbNo_dB,berTeoricoQam16Ray,'-O');
hold on;
semilogy(EbNo_dB,BerSimulado,'r-*');
grid on;