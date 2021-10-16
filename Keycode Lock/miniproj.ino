#include <Keypad.h>

 const byte rows = 4; //four rows
 const byte cols = 4; //four columns
 char keys[rows][cols] = {
 {'1','2','3','A'},
 {'4','5','6','B'},
 {'7','8','9','C'},
 {'*','0','#','D'}
 };
 byte rowPins[rows] = {0,1,2,3}; //connect to the row pinouts of the keypad
 byte colPins[cols] = {4,5,6,7}; //connect to the column pinouts of the keypad
 
  
 Keypad keypad = Keypad( makeKeymap(keys), rowPins, colPins, rows, cols );
  
 char keycode[4] = {'1','2','3','4'};
 char entry[4] = {'0','0','0','0'};
 int count = 0;
 unsigned long timer = 20000;
 bool timeout;
 bool correct = false;
  
void setup() {
  // put your setup code here, to run once
  pinMode(12, OUTPUT);
  pinMode(13, OUTPUT);

  for (int i = 0; i < 8; i++)
  {
    pinMode(i, INPUT);
  }
  keypad.setHoldTime(50);
  keypad.setDebounceTime(50);
}

void loop() {
  // put your main code here, to run repeatedly:
  if(!correct){
    digitalWrite(12, HIGH);
    digitalWrite(13, LOW);
  }
  timeout = true;
  
  for (int i = 0; i < 3; i++)
  {
    entry[i] = entry[i+1];
  }
  
  while (millis() < timer)
  {
    entry[3] = keypad.getKey();
    
    if (entry[3] != NULL)
    {
      count++;
      correct = false;
      timeout = false;
      break;
    }
  }

  if (timeout){
    count = 0;
  }
  else {
    if (count == 4){
      if ((entry[0] == keycode[0]) && (entry[1] == keycode[1]) && (entry[2] == keycode[2]) && (entry[3] == keycode[3])){
        digitalWrite(13, HIGH);
        digitalWrite(12, LOW);
        correct = true;
      }
      else{
        count = 0;
      }
    }
  }
  timer = millis() + 10000;
}
/*
#include <Keypad.h>

 const byte rows = 4; //four rows
 const byte cols = 4; //four columns
 char keys[rows][cols] = {
 {'1','2','3','A'},
 {'4','5','6','B'},
 {'7','8','9','C'},
 {'*','0','#','D'}
 };
 byte rowPins[rows] = {0,1,2,3}; //connect to the row pinouts of the keypad
 byte colPins[cols] = {4,5,6,7}; //connect to the column pinouts of the keypad
 
 Keypad keypad = Keypad( makeKeymap(keys), rowPins, colPins, rows, cols );
 
 char keycode[4] = {'1','2','3','4'};
 char entry[4] = {'0','0','0','0'};
 int count = 1;
 
void setup() {
  // put your setup code here, to run once
  pinMode(12, OUTPUT);
  pinMode(13, OUTPUT);

  for (int i = 0; i < 8; i++)
  {
    pinMode(i, INPUT);
  }
  keypad.setHoldTime(50);
  keypad.setDebounceTime(50);
  //Serial.begin(9600);
}

void loop() {
  // put your main code here, to run repeatedly:
  entry[0] = entry[1];
  entry[1] = entry[2];
  entry[2] = entry[3];
 
  entry[3] = keypad.waitForKey();
  //Serial.print(entry[3]);
 
  if (count == 4){
    if ((entry[0] == keycode[0]) && (entry[1] == keycode[1]) && (entry[2] == keycode[2]) && (entry[3] == keycode[3])){
      digitalWrite(12, HIGH);
      delay(2000);
      digitalWrite(12, LOW);
    }
    else {
      digitalWrite(13, HIGH);
      delay(2000);
      digitalWrite(13, LOW);
    }
    count = 1;
  }
  else{
    count++;
  }
}*/
