clear all;
close all;

% =========================================
% Partie I : Etude de la chaine DVBS-S    =
% dans un canal Gaussien                  =
%                                         =
%                                         =
% 2 Codage canal                          =
%                                         =
% 2.3 Ajout du code de Reed Salomon       =
% SANS ENTRELACEUR                        =
%                                         =
% =========================================

% Constantes
Fe = 24000;
Rb = 3000;
Tb = 1/Rb;
SPAN = 4;

% Modulation QPSK - 4 symboles
n = 2;    % nb de bits par symboles
M = 2^n;  % span 4 symbols car QPSK

Rs = Rb/n;
Ns = Fe/Rs; % nbe échantillons / symboles
Ts = 1/Rs;
Te = 1/Fe;

% Param codage RS
Nrs = 204;
Krs = 188;

% Génération aléatoire des N bits = information
N = Krs*Nrs*2;
bits = randi([0 1], 1, N);

%% Ajout du code de Reed Salomon %%%%%%%%%%%%%%%%%%%%

% Ajout du codage externe de Reed-Solomon
rs_encoder = comm.RSEncoder(204,188, 'BitInput',true);
bits_RS = step(rs_encoder, bits.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Code convoutif (3,1/2)
% Définition des polynomes générateurs
g1 = 5; % 5 octal -> 101
g2 = 7; % 7 octal -> 111

% Poinconnage avec la matrice P
P = [1 1 0 1];

% Création du treillis pour ce code convolutif (3,1/2)
trellis = poly2trellis(3, [g1 g2]); %Expliquer
bits_code = convenc(bits_RS, trellis, P);  %Expliquer 


% Mapping des bits en symboles QPSK
% Génération des symboles
% On passe les bits en +1 -1
symboles_bits = bits_code*2-1;

%  1 + j = 00 = -1 -1
% -1 + j = 10 =  1 -1
% -1 - j = 11 =  1  1
%  1 - j = 01 = -1  1

% Génération des symboles
symboles = zeros(1, length(bits_code)/2);
for i = 1:2:length(bits_code)
    symboles((i+1)/2) = -symboles_bits(i) - 1i * symboles_bits(1+i);
end


% Emetteur
signal = kron(symboles, [1 zeros(1, Ns -1)]);   % On prolonge de symbole avec N-1 zéros derrière

% Filtre de mise en forme 
% Racine de cosinus surélevé Beta = 0.35
Beta = 0.35;
h_n = rcosdesign(Beta, SPAN, Ns,"sqrt"); 

% Signal en sortie de l'émetteur
x = conv(signal, h_n, 'full');

% % Tracé du signal modulé
% figure;
% plot(real(x), 'b');
% hold on
% plot(imag(x), 'r');
% title('Signal modulé (partie réelle et imaginaire)');
% xlabel('Temps');
% ylabel('Amplitude');
% legend('Partie réelle', 'Partie imaginaire');
% grid on;
% hold off;

% ajout d'un canal AWGN à la chaine de transmission
Eb_sur_N0 = -8:4;   % Eb/N0 variant de 0dB à 4dB pas pas de 1dB

% Je ne pense pas que la boucle for soit encore correcte car quand on fait
% varier le SNR de 0 à 4 on ne il n'y a pas de décroissance du TEB et TES

% Création de la liste TEB pour les comparer
mat_TEB_hard = zeros(1,length(Eb_sur_N0));
mat_TEB_soft = zeros(1,length(Eb_sur_N0));

for k = 1:length(Eb_sur_N0)

    % Calcul de la puissance du bruit en fonction du SNR
    Eb_sur_N0_lineaire= 10^(Eb_sur_N0(k) / 10); % Conversion du SNR en échelle linéaire
    
    % Calculez la puissance du signal à partir du signal modulé (x)
    variance_bruit = (var(x)*Ns)/(2*log2(M)*Eb_sur_N0_lineaire);
    
    % Générez un bruit gaussien pour les parties réelles et imaginaires
    bruit_reel = sqrt(variance_bruit / 2) * randn(size(x));
    bruit_imag = sqrt(variance_bruit / 2) * randn(size(x));

    % signal à l'entrée du récepteur
    x_bruite = x + bruit_reel + 1i*bruit_imag;


    % Récepteur

    % avec le filtre adapté
    % donc on inverse la réponse impulsionnelle du filtre en émission
    z = conv(x_bruite, fliplr(h_n), 'full');  % on inverse le filtre de mise en forme pour optenir le filtre adapté


    % % Tracé de la constellation en réception
    % % On prend tous les Ns/2 symboles
    % figure;
    % scatter(real(z(128/2 + 1:Ns:end-128/2)), imag(z(128/2 + 1:Ns:end-128/2)), 'x');
    % % scatter(real(z(:)), imag(z(:)), 'x');
    % title(sprintf('Constellation QPSK (Après transmission) EB/N0 = %f', Eb_sur_N0(k)));
    % xlabel('Partie réelle');
    % ylabel('Partie imaginaire');
    % grid on;
    
    % on enlève les 128 échantillons de retard ajouté par notre conv
    z = z(128/2 + 1:end - 128/2);
    
    % % Tracer du diagramme de l'oeil pour connaitre le moment optimal de la
    % % décision
    % eyediagram(z(1:end),Ns*2); % On regarde 2 symboles par fenetre
    % title(sprintf("Diagramme de l'oeil pour %f", Eb_sur_N0(k)));

    % échantillonnage au rythme du symbole
    z = z(:,1:Ns:end);
     
   
    % Décision sur les parties réelles et imaginaires
    % Reconstruction des symboles estimés à partir des décisions R et I
    symboles_estimes = ((real(z) > 0)*2-1) + 1i *((imag(z) > 0)*2-1); 
   

    % Décision sur la partie imaginaire et réel direct
    symbole_recu = zeros(1,2*length(symboles_estimes));
    for i = 1:2:2*length(symboles_estimes)
        symbole_recu(i) = -real(symboles_estimes((i+1)/2));
        symbole_recu(i+1) = -imag(symboles_estimes((i+1)/2));
    end
    
    bits_recu = (symbole_recu + 1) /2;

    
    % Hard
    % On se peut donc maitenant reconstituer l'information envoyé
    % Décodage du code convolutif avec l'algorithme de Viterbi Hard
    bits_decode_hard = vitdec(bits_recu, trellis, 5, 'trunc', 'hard',P);


    % Soft
    % Décodage du code convolutif avec l'algorithme de Viterbi Soft
    % On crée notre vecteur alterné reel - imag
    Vect_soft = zeros(1, 2*length(z));
    Vect_soft(1:2:end) = real(z);
    Vect_soft(2:2:end) = imag(z);

    % Décodage Viterbi Soft
    bits_decode_soft = vitdec(Vect_soft, trellis, 5, 'trunc', 'unquant', P);



    % Decodage Reed-Salomon RS(204,188)
    rs_decoder = comm.RSDecoder(Nrs,Krs);

    % HARD
    bits_decode_RS_hard = step(rs_decoder, bits_decode_hard.');

    % SOFT
    bits_decode_RS_soft = step(rs_decoder, bits_decode_soft.');
    
    % Calcul TEB 
    % Ajoute de TEB dans la matrice
    mat_TEB_hard(k) = sum(bits_decode_RS_hard.' ~= bits) / N;
    mat_TEB_soft(k) = sum(bits_decode_RS_soft.' ~= bits) / N;
    
end

% Tracer TES et TEB sur le meme graphique
figure;
plot(Eb_sur_N0, mat_TEB_soft, 'b', 'DisplayName', 'TEB_Soft')
hold on;
plot(Eb_sur_N0, mat_TEB_hard, 'r', 'DisplayName', 'TEB_Hard');
hold off;

title("Taux d'erreur en fonction du SNR");
xlabel('SNR (dB)');
ylabel('Taux d''erreur');
legend('TEB Soft','TEB Hard');
grid on;