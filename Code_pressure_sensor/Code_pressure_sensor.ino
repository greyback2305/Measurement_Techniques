#include <SPI.h>

//Chip select pin
 const int CS_PIN = 10;
 //Sensor read command 
const int SENSOR_READ_CMD =0x3F;


void setup(){

  Serial.begin(9600);

  pinMode(CS_PIN,OUTPUT);   

  SPI.begin();

  SPI.beginTransaction(SPISettings(800000, MSBFIRST, SPI_MODE0));

} 


void loop(){ 

  //Select the sensor 
  digitalWrite(CS_PIN,LOW); 

  // Send the read command
  //SPI.transfer(SENSOR_READ_CMD); 
  delayMicroseconds(10);

  //Wait for conversion 
  byte sensor_data[3];
    for (int i = 0; i < 3 ;i++) { 
      //Read the data 
      sensor_data[i]= SPI.transfer(0);
  } 

  //Deselect the sensor
  digitalWrite(CS_PIN,HIGH);

  // Convert the raw sensor data to pressure 
  int pressure_counts = ((sensor_data[0] & 0x3F) << 8) | sensor_data[1];
  //Serial.println(pressure_counts);
  float pressure = 6894.76*(((pressure_counts - 1638.00)/ (14745.00-1638.00))*2-1.00);
  //Serial.print("Pressure:"); 
  Serial.print(pressure);
  Serial.println("Pa");

  delay(100);
  //Wait for 1 second before the next reading 
}

// Coef de convertion k a calibrer : k = pression_connue / (valeur_mesuree); (pression_connue calculée theoriquement)

