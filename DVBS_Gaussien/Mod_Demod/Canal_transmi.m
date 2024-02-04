clear all;
close all;

% =========================================
% Partie I : Etude de la chaine DVBS-S    =
% dans un canal Gaussien                  =
%                                         =
%                                         =
% 1 Modulateur/ Démodulateur              =
%                                         =
% 1.3 Canal de transmission               =
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

% Génération aléatoire des n bits = information
nb_bits = 2000;
bits = randi([0 1], 1, nb_bits);

% Mapping des bits en symboles QPSK
% Génération des symboles
% On passe les bits en +1 -1
bits =bits*2-1;

%  1 + j = 00 = -1 -1
% -1 + j = 10 =  1 -1
% -1 - j = 11 =  1  1
%  1 - j = 01 = -1  1

% Génération des symboles
symboles = zeros(1, nb_bits/2);
for i = 1:2:nb_bits
    symboles((i+1)/2) = -bits(i) - 1i * bits(1+i);
end

% % Tracé de la constellation de l'émetteur
% figure;
% scatter(real(symboles), imag(symboles), 'filled');
% title('Constellation QPSK (Avant transmission)');
% xlabel('Partie réelle');
% ylabel('Partie imaginaire');
% grid on;

% Emetteur
signal = kron(symboles, [1 zeros(1, Ns -1)]);   % On prolonge de symbole avec N-1 zéros derrière

% Filtre de mise en forme 
% Racine de cosinus surélevé Beta = 0.35
Beta = 0.35;
h_n = rcosdesign(Beta, SPAN, Ns,"sqrt"); 

% Signal en sortie de l'émetteur
x = conv(signal, h_n, 'full');


%% Canal de transmission
% Ajout d'un canal AWGN à la chaine de transmission
Eb_sur_N0 = 0:4;   % Eb/N0 variant de 0dB à 4dB pas pas de 1dB


% Création des matrices TEB et TES pour les comparer
mat_TEB = zeros(1,length(Eb_sur_N0));
mat_TES = zeros(1,length(Eb_sur_N0));

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


    % Tracé de la constellation en réception
    % On prend tous les Ns/2 symboles
    figure;
    scatter(real(z(128/2 + 1:Ns:end-128/2)), imag(z(128/2 + 1:Ns:end-128/2)), 'x');
    % scatter(real(z(:)), imag(z(:)), 'x');
    title(sprintf('Constellation QPSK (Après transmission) EB/N0 = %f', Eb_sur_N0(k)));
    xlabel('Partie réelle');
    ylabel('Partie imaginaire');
    grid on;
    
    % on enlève les 128 échantillons de retard
    z = z(128/2 + 1:end - 128/2);
    
    % Tracer du diagramme de l'oeil pour connaitre le moment optimal de la
    % décision
    % eyediagram(z(Ns*2:end-Ns*2),Ns*2); % On regarde 2 symboles par fenetre
    % title(sprintf("Diagramme de l'oeil pour %f", Eb_sur_N0(k)));
    
    % Décision sur les parties réelles et imaginaires
    decision_R = (real(z(:,1:Ns:end)) < 0); % Décision sur la partie réelle
    decision_I = (imag(z(:,1:Ns:end)) < 0); % Décision sur la partie imaginaire
    

    
    % Décision sur les parties réelles et imaginaires
    % Reconstruction des symboles estimés à partir des décisions R et I
    symboles_estimes = ((real(z(:,1:Ns:end)) > 0)*2-1) + 1i *((imag(z(:,1:Ns:end)) > 0)*2-1); 
    % Décision sur la partie imaginaire et réel direct
    
    % On se peut donc maitenant reconstituer l'information envoyé
    bits_recu = zeros(1, nb_bits); % Initialisation des bits reçus
    
    for i = 1:2:nb_bits
        bits_recu(i) = -real(symboles_estimes((i+1)/2));
        bits_recu(i+1) = -imag(symboles_estimes((i+1)/2));
    end
    
    
    % Calcul du taux d'erreur de symbole (TES)
    % Ajout de TES dans la matrice 
    mat_TES(k) = sum(symboles_estimes ~= symboles) / length(symboles);
    
    % Calcul TEB 
    % Ajoute de TEB dans la matrice
    mat_TEB(k) = sum(bits_recu ~= bits) / nb_bits;
    
end

% Tracer TES et TEB sur le meme graphique
figure;
hold on;
plot(Eb_sur_N0, mat_TES, 'b', 'DisplayName', 'TES');
plot(Eb_sur_N0, mat_TEB, 'r', 'DisplayName', 'TEB');
hold off;

title("Taux d'erreur en fonction du SNR");
xlabel('SNR (dB)');
ylabel('Taux d''erreur');
legend('TES', 'TEB');
grid on;