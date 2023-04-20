

# # # # # # # # # #
# # Beginn des Scenario Files
# # # # # # # # # #


# # # # # # # # # #
# # Header
# # # # # # # # # #


# Definition der beteiligten Dateien
scenario = "TempInt_Demask_s3.scl";
pcl_file = "TempInt_Demask_s3.pcl";

# Definition der Button-Code-Parameter
response_matching = simple_matching;
active_buttons = 6;
button_codes = 91, 92, 93, 94, 95, 96;
default_background_color = 255, 255, 255; 


# Definition der Port-Code-Parameter
write_codes = true;
pulse_width = 1;


# # # # # # # # # #
# # Beginn des Trial-Baukastens


begin;


# # # # # # # # # #
# # Definition möglicher Stimuli

# Bild des Fixations-Punktes
bitmap{
	
	filename = "fixation_point.png";
	transparent_color = 255, 255, 255;
	
} BITMAP_FIXATION_POINT;

# Bild des Target-Quadrats
bitmap{
	
	filename = "target_square.png";
	transparent_color = 255, 255, 255;
	
} BITMAP_TARGET_SQUARE;

# Bild der 1 px Quadrat-Maske
bitmap{ 
	
	filename = "mask_square_1px.png";
	transparent_color = 255, 255, 255;

} BITMAP_MASK_1;

# Bild der 2 px Quadrat-Maske
bitmap{ 
	
	filename = "mask_square_2px.png";
	transparent_color = 255, 255, 255;

} BITMAP_MASK_2;

# Bild der 3 px Quadrat-Maske
bitmap{ 
	
	filename = "mask_square_3px.png";
	transparent_color = 255, 255, 255;

} BITMAP_MASK_3;

# Bild der 4 px Quadrat-Maske
bitmap{ 
	
	filename = "mask_square_4px.png";
	transparent_color = 255, 255, 255;

} BITMAP_MASK_4;

# Bild der 5 px Quadrat-Maske
bitmap{ 
	
	filename = "mask_square_5px.png";
	transparent_color = 255, 255, 255;

} BITMAP_MASK_5;

# Bild der 6 px Quadrat-Maske
bitmap{ 
	
	filename = "mask_square_6px.png";
	transparent_color = 255, 255, 255;

} BITMAP_MASK_6;

# Bild der 7 px Quadrat-Maske
bitmap{ 
	
	filename = "mask_square_7px.png";
	transparent_color = 255, 255, 255;

} BITMAP_MASK_7;

# Bild der 8 px Quadrat-Maske
bitmap{ 
	
	filename = "mask_square_8px.png";
	transparent_color = 255, 255, 255;

} BITMAP_MASK_8;

# Bild der 9 px Quadrat-Maske
bitmap{ 
	
	filename = "mask_square_9px.png";
	transparent_color = 255, 255, 255;

} BITMAP_MASK_9;

# Bild der 10 px Quadrat-Maske
bitmap{ 
	
	filename = "mask_square_10px.png";
	transparent_color = 255, 255, 255;

} BITMAP_MASK_10;

# Bild der 11 px Quadrat-Maske
bitmap{ 
	
	filename = "mask_square_11px.png";
	transparent_color = 255, 255, 255;

} BITMAP_MASK_11;

# Bild der 12 px Quadrat-Maske
bitmap{ 
	
	filename = "mask_square_12px.png";
	transparent_color = 255, 255, 255;

} BITMAP_MASK_12;

# Bild der 13 px Quadrat-Maske
bitmap{ 
	
	filename = "mask_square_13px.png";
	transparent_color = 255, 255, 255;

} BITMAP_MASK_13;

# Bild der 14 px Quadrat-Maske
bitmap{ 
	
	filename = "mask_square_14px.png";
	transparent_color = 255, 255, 255;

} BITMAP_MASK_14;

# Bild der 15 px Quadrat-Maske
bitmap{ 
	
	filename = "mask_square_15px.png";
	transparent_color = 255, 255, 255;

} BITMAP_MASK_15;

# Bild der 16 px Quadrat-Maske
bitmap{ 
	
	filename = "mask_square_16px.png";
	transparent_color = 255, 255, 255;

} BITMAP_MASK_16;

# Bild der Antwort-Fixation
bitmap{ 
	
	filename = "fixation_response_window.png";
	transparent_color = 255, 255, 255;

} BITMAP_RESPONSE;


# # # # # # # # # #
# # Text-Defintionen

# Text der Willkommens-Folie
text{

	caption = "Liebe*r Versuchsteilnehmer*in,
	
	willkommen zurück und nochmals vielen Dank für deine Bereitschaft an diesem Experiment teilzunehmen!
	
	Viel Spaß!
	
	Drücke die LEERTASTE um zur nächsten Folie zu gelangen.";
   font ="Arial";
   font_size = 12;
   font_color = 0, 0, 0;
	max_text_width = 500;
	formatted_text = true;

} TEXT_WILLKOMMEN;

# Text zur Rekapitulation der Aufgabe
text{

	caption = "Heute werden dir wieder quadratische Zielreize und Masken in schneller Abfolge präsentiert. Deine Aufgabe ist es die Sichtbarkeit des Zielreizes (sichtbar vs. nicht sichtbar) und die zeitliche Abfolge von Zielreiz und Maske (gleichzeitig vs. nacheinander) einzuschätzen.
	
	Deinen Wahrnehmungseindruck kannst du stets am Ende eines Durchgangs abgeben sobald sich der FixationsPUNKT in ein FixationsKREUZ verwandelt hat. Führe deinen Tastendruck bitte nicht vorher aus, da sonst die Hirnströme mit Störsignalen belastet werden!

	Drücke die LEERTASTE um zur nächsten Folie zu gelangen.";
	
   font ="Arial";
   font_size = 12;
   font_color = 0, 0, 0;
	max_text_width = 500;
	formatted_text = true;

} TEXT_TASK_INFO;

# Text zur Einläutung der Hauptblöcke
text{

	caption = "Im Folgenden werden dir 8 Blöcke mit je 96 Durchgängen präsentiert. Zwischen jedem Block steht es dir frei eine Pause zu machen und beispielsweise etwas zu trinken oder die Augen zu entspannen. Versuche dies allerdings bitte auf die Pausen zu beschränken.
	
	Hast du Fragen? Dann wende dich jetzt an deine*n Experimentator*in!

	Durch Drücken der LEERTASTE startet das eigentliche Experiment.";
	
   font ="Arial";
   font_size = 12;
   font_color = 0, 0, 0;
	max_text_width = 500;
	formatted_text = true;

} TEXT_BLOCK_INFO;

# Text für Vorbereitung der Fixation
text{

	caption = "Fokussiere nun bitte den Fixationspunkt und verlasse diesen nur, wenn du deinen Wahrnehmungseindruck abgibst oder du dich in einer Pause zwischen Blöcken befinden.

	Durch Drücken der LEERTASTE signalisierst du deine Bereitschaft Block X/8 zu starten. Der Block beginnt dann sobald dein EEG frei von Störsignalen ist.]";
	
	font = "Arial";
	font_size = 12;
	font_color = 0, 0, 0; 
	max_text_width = 500;
	formatted_text = true;

} TEXT_FIXATION_PREPARATION;

# Text zur Beruhigung des EEG
text{

	caption = " ";
	
	font = "Arial";
	font_size = 12;
	font_color = 0, 0, 0; 
	max_text_width = 500;
	formatted_text = true;

} TEXT_EEG_PREPARATION;

# Text für die Pausen zwischen den Blöcken
text{

	caption = "Du hast Block X/8 erfolgreich abgeschlossen!
	
	Durch Drücken der LEERTASTE startest du den nächsten Block.";
	
	font = "Arial";
	font_size = 12;
	font_color = 0, 0, 0; 
	max_text_width = 500;
	formatted_text = true;

} TEXT_PAUSE;

# Text zum Beenden des Experiments
text{
	
	caption = "Du hast diese Sitzung erfolgreich abgeschlossen!
	
	Vielen Dank für deine Teilnahme!
	
	Signalisiere bitte deinem*r Experimentator*in durch Drücken der LEERTASTE, dass du das Experiment beendet hast. Du wirst dann in Kürze aus der Dunkelheit befreit.";
	
	font = "Arial";
	font_size = 12;
	font_color = 0, 0, 0;
	max_text_width = 500;
	formatted_text = true;

} TEXT_ENDE;


# # # # # # # # # #
# # Picture-Defintionen

# Definition der Willkommens-Folie
picture{

	text TEXT_WILLKOMMEN;
	x = 0;
	y = 0;

} PICTURE_WILLKOMMEN;

# Definition der Aufgaben-Instruktions-Folie
picture{

	text TEXT_TASK_INFO;
	x = 0;
	y = 0;

} PICTURE_TASK_INFO;

# Definition der Instruktions-Folie zum Starten der Blöcke
picture{

	text TEXT_BLOCK_INFO;
	x = 0;
	y = 0;

} PICTURE_BLOCK_INFO;

# Definition der Vorbereitung auf die Fixation Folie
picture{

	text TEXT_FIXATION_PREPARATION;
	x = 0;
	y = -125;
	
	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;
	
} PICTURE_FIXATION_PREPARATION;

# Definition der EEG-Beruhigungs-Folie
picture{

	text TEXT_EEG_PREPARATION;
	x = 0;
	y = -125;
	
	bitmap BITMAP_RESPONSE;
	x = 0;
	y = 0;
	
} PICTURE_EEG_PREPARATION;

# Definition des System-Pausieren-Bildschirms zum Starten des EEG
picture{

	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;

} PICTURE_SYSTEM_BREAK;

# Definition der Präsentation des Fixationspunktes
picture{

	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;
	
} PICTURE_FIXATION;

# Definition der quadratischen Target-Präsentation
picture{

	bitmap BITMAP_TARGET_SQUARE;
	x = 0;
	y = 0;
	
} PICTURE_TARGET_SQUARE;

# Definition der 1px quadratischen Masken-Präsentation
picture{

	bitmap BITMAP_MASK_1;
	x = 0;
	y = 0;
	
	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;

} PICTURE_MASK_1;

# Definition der 2px quadratischen Masken-Präsentation
picture{

	bitmap BITMAP_MASK_2;
	x = 0;
	y = 0;
	
	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;

} PICTURE_MASK_2;

# Definition der 3px quadratischen Masken-Präsentation
picture{

	bitmap BITMAP_MASK_3;
	x = 0;
	y = 0;
	
	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;

} PICTURE_MASK_3;

# Definition der 4px quadratischen Masken-Präsentation
picture{

	bitmap BITMAP_MASK_4;
	x = 0;
	y = 0;
	
	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;

} PICTURE_MASK_4;

# Definition der 5px quadratischen Masken-Präsentation
picture{

	bitmap BITMAP_MASK_5;
	x = 0;
	y = 0;
	
	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;

} PICTURE_MASK_5;

# Definition der 6px quadratischen Masken-Präsentation
picture{

	bitmap BITMAP_MASK_6;
	x = 0;
	y = 0;
	
	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;

} PICTURE_MASK_6;

# Definition der 7px quadratischen Masken-Präsentation
picture{

	bitmap BITMAP_MASK_7;
	x = 0;
	y = 0;
	
	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;

} PICTURE_MASK_7;

# Definition der 8px quadratischen Masken-Präsentation
picture{

	bitmap BITMAP_MASK_8;
	x = 0;
	y = 0;
	
	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;

} PICTURE_MASK_8;

# Definition der 9px quadratischen Masken-Präsentation
picture{

	bitmap BITMAP_MASK_9;
	x = 0;
	y = 0;
	
	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;

} PICTURE_MASK_9;

# Definition der 10px quadratischen Masken-Präsentation
picture{

	bitmap BITMAP_MASK_10;
	x = 0;
	y = 0;
	
	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;

} PICTURE_MASK_10;

# Definition der 11px quadratischen Masken-Präsentation
picture{

	bitmap BITMAP_MASK_11;
	x = 0;
	y = 0;
	
	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;

} PICTURE_MASK_11;

# Definition der 12px quadratischen Masken-Präsentation
picture{

	bitmap BITMAP_MASK_12;
	x = 0;
	y = 0;
	
	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;

} PICTURE_MASK_12;

# Definition der 13px quadratischen Masken-Präsentation
picture{

	bitmap BITMAP_MASK_13;
	x = 0;
	y = 0;
	
	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;

} PICTURE_MASK_13;

# Definition der 14px quadratischen Masken-Präsentation
picture{

	bitmap BITMAP_MASK_14;
	x = 0;
	y = 0;
	
	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;

} PICTURE_MASK_14;

# Definition der 15px quadratischen Masken-Präsentation
picture{

	bitmap BITMAP_MASK_15;
	x = 0;
	y = 0;
	
	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;

} PICTURE_MASK_15;

# Definition der 16px quadratischen Masken-Präsentation
picture{

	bitmap BITMAP_MASK_16;
	x = 0;
	y = 0;
	
	bitmap BITMAP_FIXATION_POINT;
	x = 0;
	y = 0;

} PICTURE_MASK_16;


# Definition des Antwort-Fensters
picture{

	bitmap BITMAP_RESPONSE;
	x = 0;
	y = 0;

} PICTURE_RESPONSE;

# Definition der Pausen-Folie zwischen den Durchgängen
picture{

	text TEXT_PAUSE;
	x = 0;
	y = 0;

} PICTURE_PAUSE;

# Definition der Folie zum Abschluss des Experiments
picture{
	
	text TEXT_ENDE;
	
	x = 0;
	y = 0;

} PICTURE_ENDE;

# Definition der Präsentation des Fixationspunktes als Default-Image
picture{

	bitmap BITMAP_FIXATION_POINT;
	
	x = 0;
	y = 0;
	
} default;


# # # # # # # # # #
# # Trial-Definitionen

# Definition des Willkommens-Trials
trial{

	trial_duration = forever;
	trial_type = specific_response;
	terminator_button = 5;

	picture PICTURE_WILLKOMMEN;
	
} TRIAL_WILLKOMMEN;

# Definition des Trials zur Aufgaben-Instruktion
trial{

	trial_duration = forever;
	trial_type = specific_response;
	terminator_button = 5;

	picture PICTURE_TASK_INFO;
	
} TRIAL_TASK_INFO;

# Definition des Trials kurz vor dem Start der Hauptblöcke
trial{

	trial_duration = forever;
	trial_type = specific_response;
	terminator_button = 5;

	picture PICTURE_BLOCK_INFO;
	
} TRIAL_BLOCK_INFO;

# Definition des Trials zur Vorbereitung der Fixation
trial{

	trial_duration = forever;
	trial_type = specific_response;
	terminator_button = 5;

	picture PICTURE_FIXATION_PREPARATION;
	
} TRIAL_FIXATION_PREPARATION;

# Definition des Trials zur Beruhigung der EEG-Kurven
trial{

	trial_duration = forever;
	trial_type = specific_response;
	terminator_button = 6;

	picture PICTURE_EEG_PREPARATION;
	
} TRIAL_EEG_PREPARATION;

# Definition des Trials zur Pausenpräsentation zwischen den Blöcken
trial{

	trial_duration = forever;
	trial_type = specific_response;
   terminator_button = 5;

	stimulus_event {
		
		picture PICTURE_PAUSE;       
		
		code = "PAUSE";
		
	} STIMULUS_EVENT_PAUSE;
	
} TRIAL_PAUSE;

# Definition eines Trials, welches dem PC vor und nach jedem Block Zeit gibt um die Aufnahme des EEG zu starten bzw. zu beenden
trial{

	trial_duration = 997;
	all_responses = false;
	
	picture PICTURE_SYSTEM_BREAK;
	
} TRIAL_SYSTEM_BREAK;

# Definition des eigentlichen Trials
trial{
	
	trial_duration = stimuli_length;
	all_responses = false;
	
	stimulus_event{
		
		picture PICTURE_FIXATION;
		
		response_active = false;
		duration = next_picture;
		port_code = 10;
		code = "FIXATION_FIRST";		
		
	} STIMULUS_EVENT_FIXATION_FIRST;
	
	stimulus_event{
		
		picture PICTURE_TARGET_SQUARE;
		
		response_active = false;
		duration = 22;
		time = 997;
		port_code = 108;
		code = "TARGET";		
		
	} STIMULUS_EVENT_TARGET;
	
	stimulus_event{
		
		picture PICTURE_MASK_8; 
		
		response_active = false;
		duration = 30;
		deltat = 22;
		port_code = 208;
		code = "MASK";		
		
	} STIMULUS_EVENT_MASK;
	
	stimulus_event{
		
		picture PICTURE_FIXATION;
		
		response_active = false;
		duration = 1003;
		deltat = 30;
		code = "FIXATION_SECOND";
		
	} STIMULUS_EVENT_FIXATION_SECOND;
	
} TRIAL_MAIN;

# Definition des Trials zur Antwortabgabe
trial{

	trial_duration = forever;
	trial_type = specific_response;
	terminator_button = 1,2,3,4;

	stimulus_event{

		picture PICTURE_RESPONSE;
		
		response_active = true;
		code = "RESPONSE";
		port_code = 90;
		
	} STIMULUS_EVENT_RESPONSE;

} TRIAL_RESPONSE;

# Definition des Trials zum Abschluss des Experiments
trial{

	trial_duration = forever;
	trial_type = specific_response;
	terminator_button = 5;

	picture PICTURE_ENDE;
	
} TRIAL_ENDE;


# # # # # # # # # #
# # Array-Definitionen

# Definition des Arrays, aus welchem die Masken ausgewählt werden
array{
	
	picture PICTURE_MASK_1;
	
	picture PICTURE_MASK_2;
	
	picture PICTURE_MASK_3;
	
	picture PICTURE_MASK_4;
	
	picture PICTURE_MASK_5;
	
	picture PICTURE_MASK_6;
	
	picture PICTURE_MASK_7;
	
	picture PICTURE_MASK_8;
	
	picture PICTURE_MASK_9;
	
	picture PICTURE_MASK_10;
	
	picture PICTURE_MASK_11;
	
	picture PICTURE_MASK_12;
	
	picture PICTURE_MASK_13;
	
	picture PICTURE_MASK_14;
	
	picture PICTURE_MASK_15;
	
	picture PICTURE_MASK_16;

} ARRAY_MASKS;


# # # # # # # # # #
# # Ende des Scenario Files
# # # # # # # # # #

