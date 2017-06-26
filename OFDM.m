function varargout = OFDM(varargin)
% OFDM MATLAB code for OFDM.fig
%      OFDM, by itself, creates a new OFDM or raises the existing
%      singleton*.
%
%      H = OFDM returns the handle to a new OFDM or the handle to
%      the existing singleton*.
%
%      OFDM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OFDM.M with the given input arguments.
%
%      OFDM('Property','Value',...) creates a new OFDM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OFDM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OFDM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OFDM

% Last Modified by GUIDE v2.5 21-Jun-2017 12:49:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OFDM_OpeningFcn, ...
                   'gui_OutputFcn',  @OFDM_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before OFDM is made visible.
function OFDM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to OFDM (see VARARGIN)

% Choose default command line output for OFDM
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%Inicialização de Variáveis
global IntGuard; %CheckBox de Intervalor de Guarda
IntGuard=0; %Inicia Off
global Ocp; %Opção Prefixo Cíclico
Ocp=0;
global Ogi; %Opção Intervalo de Guarda
Ogi=0;
global Pcp; %Porcentagem Intervalo de Guarda
Pcp=25;

% UIWAIT makes OFDM wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = OFDM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function EbNoIn_Callback(hObject, eventdata, handles)
% hObject    handle to EbNoIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EbNoIn as text
%        str2double(get(hObject,'String')) returns contents of EbNoIn as a double


% --- Executes during object creation, after setting all properties.
function EbNoIn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EbNoIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function EbNoFim_Callback(hObject, eventdata, handles)
% hObject    handle to EbNoFim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EbNoFim as text
%        str2double(get(hObject,'String')) returns contents of EbNoFim as a double


% --- Executes during object creation, after setting all properties.
function EbNoFim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EbNoFim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Qsp_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Qsp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Qsp_edit as text
%        str2double(get(hObject,'String')) returns contents of Qsp_edit as a double


% --- Executes during object creation, after setting all properties.
function Qsp_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Qsp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Qinf_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Qinf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Qinf_edit as text
%        str2double(get(hObject,'String')) returns contents of Qinf_edit as a double


% --- Executes during object creation, after setting all properties.
function Qinf_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Qinf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in chBox.
function chBox_Callback(hObject, eventdata, handles)
% hObject    handle to chBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns chBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from chBox


% --- Executes during object creation, after setting all properties.
function chBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white'); 
end
set(hObject,'String',{'16-QAM';'BPSK'});

% --- Executes on selection change in chBox3.
function chBox3_Callback(hObject, eventdata, handles)
% hObject    handle to chBox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns chBox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from chBox3


% --- Executes during object creation, after setting all properties.
function chBox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chBox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',{'AWGN';'Rayleigh Fading'});



function Pcp_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Pcp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pcp_edit as text
%        str2double(get(hObject,'String')) returns contents of Pcp_edit as a double


% --- Executes during object creation, after setting all properties.
function Pcp_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pcp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chBox2.
function chBox2_Callback(hObject, eventdata, handles)
% hObject    handle to chBox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global IntGuard;
IntGuard=0;

% Hint: get(hObject,'Value') returns toggle state of chBox2
if (get(hObject,'Value') == get(hObject,'Max'))
	
    set(handles.BtGroup,'visible','on');
    IntGuard=1;
else
	 set(handles.BtGroup,'visible','off');
     IntGuard=0;
end

function nRep_Callback(hObject, eventdata, handles)
% hObject    handle to nRep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nRep as text
%        str2double(get(hObject,'String')) returns contents of nRep as a double


% --- Executes during object creation, after setting all properties.
function nRep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nRep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gerar.
function gerar_Callback(hObject, eventdata, handles)
st='Status: ';
set(handles.Teste,'String',st);

% Configuração e Parâmetros de entrada

NumSub=str2num(get(handles.Qsp_edit,'String'));  %Quantidade de subportadoras do sistema OFDM
Qsimb=str2num(get(handles.Qinf_edit,'String'));   %Quantidade de símbolos à ser transmitido
Nfft = 2^(nextpow2(NumSub));%Total de subportadoras usadas e não-usadas
Ncp=0; %Quantidade de prefixo ciclico ou intervalo de guarda
Nbg=Nfft-NumSub; % Banda Virtual (ZeroPading)

Qrep = str2num(get(handles.nRep,'String'));
% Qrep=100;

TipoMod=0; %Tipo de modulação
TipoCanal=0; %Tipo de canal

%Variaveis de Erros
EbNoin = str2num(get(handles.EbNoIn,'String'));
EbNofim = str2num(get(handles.EbNoFim,'String'));
EbNo_dB = EbNoin:5:EbNofim; % Variação de Eb/N0 dB


% Verificar Banda de Guarda
global IntGuard;
global Opc Ogi Pcp;

if IntGuard==1    % SE BANDA DE GUARDA ESTIVER ATIVO
    if (get(handles.rbtGI,'Value') == get(handles.rbtGI,'Max')) %Se Intervalo de guarda etiver ativo
        Ogi=1;
        Opc=0;
        set(handles.Teste,'String',[get(handles.Teste,'String') '|Intervalo de Guarda']);
        
    else  %Se Prefixo Cíclico estiver ativo
        Ogi=0;
        Opc=1;
        set(handles.Teste,'String',[get(handles.Teste,'String') '|Prefico Cíclico']);
    end
    
Pcp = str2num(get(handles.Pcp_edit,'String'));
Ncp = (Pcp/100)*Nfft;

else
    set(handles.Teste,'String',[get(handles.Teste,'String') '|Sem Int.Guarda']);
end

% Verificar Tipo de Modulação
if get(handles.chBox,'Value')==1
    TipoMod=1; %Modulação 16-QAM
    set(handles.Teste,'String',[get(handles.Teste,'String') '|16-QAM']);
elseif get(handles.chBox,'Value')==2
    TipoMod=2; %Modulação BPSK
    set(handles.Teste,'String',[get(handles.Teste,'String') '|BPSK']);
end


%Verificar Tipo de Canal
if get(handles.chBox3,'Value')==1 %Canal com Ruído AWGN
    TipoCanal=1; %AWGN
    set(handles.Teste,'String',[get(handles.Teste,'String') '|AWGN']);
elseif get(handles.chBox3,'Value')==2 %Canal com Rayleigh Fading
    TipoCanal=2; %Rayleigh
    set(handles.Teste,'String',[get(handles.Teste,'String') '|Rayleigh']);
end



%Erros

snr=zeros(1,length(EbNo_dB));
Erros=zeros(1,length(EbNo_dB));
ErrosG=zeros(1,length(EbNo_dB));


for i=1:length(EbNo_dB)
    for k=1:Qrep
        %Iniciar as Variáveis
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
    
        
%% Gerador de Bits e Modulação

        if TipoMod==1 % Modulação 16-QAM
            
            Mod = 4; %Tipo de modulação (4 bits/simbolos)
            Qbits= Mod*NumSub*Qsimb; %Quantidade de bits
            a=randi([0 1],1,Qbits); %Gerando Bits aleatórios
            ainf=bi2de(transpose(reshape(a,Mod,length(a)/Mod)),'left-msb');
        
            mod=qammod(transpose(ainf),16,0,'gray')/sqrt(10); %Modulação 16-Qam

        elseif TipoMod==2 % Modulação BPSK
            Mod = 1; %Tipo de modulação (4 bits/simbolos)
            Qbits= Mod*NumSub*Qsimb;
            mod=2*round(rand(1,Qbits))-1; %Gerando Bits em modulação BPSK
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

        vtFft=(1:Nfft); % Reiniciar variável
        
%% Adicionar Prefixo Cíclico/Intervalo de Guarda
        if IntGuard==1 && Opc==1 && Ogi==0  % Se for Prefixo Cíclico
    
            % Prefixo Cíclico
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
        Tx = xtime; % Sem Intervalo de Guarda e nem prefixo Cíclico
        
  %% Modelagem do Canal   
        
       
       snr(i)=EbNo_dB(i) + 10*log10(Mod) + 10*log10(NumSub/Nfft);
       
        if IntGuard==1 % Se tiver Intervalo de Guarda(Zero Padding) ou Prefixo Cíclico
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
        
         
         
%% Removendo Prefixo Cíclico/Intervalo de Guarda
         Rx = Rx_noise;
         
         if IntGuard==1 
             for m=1:Qsimb
                RxG = [RxG Rx_noiseG(vtRCp)];

               vtRCp=vtRCp+(Nfft+Ncp);
  
             end
         end
         
         
         %% FFT
         
         vtFft=(1:Nfft); %Reiniciando Variável
         
        %FFT sem Banda de Guarda
         for m=1:Qsimb
            yfreq = [ yfreq fft(Rx(vtFft))];
            vtFft=vtFft+Nfft;
         end
        
        
        if IntGuard==1
            %FFT com Banda de Guarda
            vtFft=(1:Nfft); %Reiniciando Variável
            
            for m=1:Qsimb
                yfreqG = [ yfreqG fft(RxG(vtFft))];
                vtFft=vtFft+Nfft;
            end
        end
        
        
        
%% Retirada da Banda de Guarda (Banda Virtual)
         vtRBg=(2:(floor(NumSub/2)+1)); %Reiniciando Variável
        vtRBgf=(floor(NumSub/2)+(Nbg+1):Nfft);
        
        for m=1:Qsimb
             ytime=[ytime [yfreq(vtRBg) yfreq(vtRBgf)]];
             vtRBg=vtRBg+Nfft;
             vtRBgf=vtRBgf+Nfft;
        end
        
        if IntGuard==1
             vtRBg=(2:(floor(NumSub/2)+1)); %Reiniciando Variável
            vtRBgf=(floor(NumSub/2)+(Nbg+1):Nfft);
            for m=1:Qsimb
                 ytimeG=[ytimeG [yfreqG(vtRBg) yfreqG(vtRBgf)]];
                 vtRBg=vtRBg+Nfft;
                 vtRBgf=vtRBgf+Nfft;
            end
        end
        
        
        
%%Equalização caso (Rayleigh)
         if TipoCanal==2 % Se tiver Espalhaento Rayleigh

            
            if IntGuard==1 %Se tiver intervalo de guarda
                
                vtRBg=(2:(floor(NumSub/2)+1)); %Reiniciando Variável
                vtRBgf=(floor(NumSub/2)+(Nbg+1):Nfft);

                hFG= fft(hG,Nfft);
                hFG=[hFG(vtRBg) hFG(vtRBgf)];
                vtEq = (1:NumSub); 

                for m = 1:Qsimb                
                ytimeG(vtEq) = ytimeG(vtEq)./hFG; %equalizador
                vtEq = vtEq + NumSub; 
                end
            end
            
            vtRBg=(2:(floor(NumSub/2)+1)); %Reiniciando Variável
            vtRBgf=(floor(NumSub/2)+(Nbg+1):Nfft);
            
            hF= fft(h,Nfft);
            hF=[hF(vtRBg) hF(vtRBgf)];
            vtEq = (1:NumSub); 
            
            for m = 1:Qsimb                
            ytime(vtEq) = ytime(vtEq)./hF; %equalizador
            vtEq = vtEq + NumSub;
            end           
         end
        
        
%% Demodulação
        
        if TipoMod==1  %Demodulação 16-QAM
            
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

        elseif TipoMod==2    %Demodulação BPSK
            
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
                errG=sum(mod~=infOutG);
                ErrosG(i)=ErrosG(i)+errG;
            end
            
            %Erros
            err=sum(mod~=infOut);
            Erros(i)=Erros(i)+err;
        end      
        
        
        
    end
end


%% Cálculo BER

EbNoLin=10.^(EbNo_dB/10); %Eb/N0 em escala linear
gamma_c=EbNoLin*Mod;

BerTeoricoBpskAwgn=(1/2)*erfc(sqrt(10.^(EbNo_dB/10))); % BER BPSK AWGN
berTeoricoBpskRay=0.5*(1-sqrt(EbNoLin./(1+EbNoLin))); %BER BPSK Rayleigh

BerTeoricoQam16Awgn = (1/Mod)*3/2*erfc(sqrt(Mod*0.1*(10.^(EbNo_dB/10))));%BER 16-QAM AWGN
berTeoricoQam16Ray = 3/8 * ( 1 - sqrt(2/5*gamma_c/Mod./(1+2/5*gamma_c/Mod)) ); %BER 16-QAM Rayleigh

BerSimulado = Erros/(Qsimb)/(Mod)/(Qrep)/(NumSub); %BER SIMULADO

if IntGuard==1 %Se tiver Int. Guarda
    BerSimuladoG = ErrosG/(Qsimb)/(Mod)/(Qrep)/(NumSub); %BER SIMULADO
end

%% Plot

figure(1);
%Plot Teórico
if TipoMod==1 %16-QAM Teórico

    semilogy(EbNo_dB,BerTeoricoQam16Awgn,'k-O'); %AWGN
    
     hold on;
    if TipoCanal==2
      semilogy(EbNo_dB,berTeoricoQam16Ray,'b-'); %Rayleigh
      hold on; 
      ylim([1*10^-7 1]);
    end
    
    
elseif TipoMod==2 %BPSK Teórico

    semilogy(EbNo_dB,BerTeoricoBpskAwgn,'-O'); %AWGN
    
     hold on;
    if TipoCanal==2
      semilogy(EbNo_dB,berTeoricoBpskRay,'b-'); %Rayleigh
      hold on;  
      ylim([1*10^-7 1]);
    end
    
    
     
end

%Plot Simulado
hold on;
semilogy(EbNo_dB,BerSimulado,'r-*'); % SIMULADO
grid on;
hold on;
if IntGuard==1 %Se tiver Int. Guarda
    semilogy(EbNo_dB,BerSimuladoG,'m->'); %SIMULADO com intervalo de guarda
end
grid on;

title('BER vs EbN0dB OFDM Teórico e Simulado');
xlabel('Eb/N0 dB');
ylabel('BER');


%% Legendas

if TipoCanal==1 
    if IntGuard==1
        legend('Teórico AWGN','Simulado AWGN S/ Int. de Guarda','Simulado AWGN C/ Int. de Guarda');
    else
         legend('Teórico AWGN','Simulado AWGN S/ Int. de Guarda');
    end
elseif TipoCanal==2
        if IntGuard==1
        legend('Teórico AWGN','Teórico Rayleigh','Simulado Rayleigh S/ Banda de Guarda','Simulado Rayleigh C/ Banda de Guarda');
        else
            legend('Teórico AWGN','Teórico Rayleigh','Simulado Rayleigh S/ Banda de Guarda');
         end
end

if TipoMod==1
    if TipoCanal==1
        st2 = texlabel(get(handles.Teste,'String')); 
        text(EbNoin+0.5,10^-2,st2);
    elseif TipoCanal==2
        st2 = texlabel(get(handles.Teste,'String')); 
        text(EbNoin+0.5,10^-3.5,st2);
    end
elseif TipoMod==2
st2 = texlabel(get(handles.Teste,'String'));  
text(EbNoin+0.5,10^(-3.5),st2);
end











    




