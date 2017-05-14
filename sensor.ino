int DHpin = 8;
byte dat [5];
byte read_data ()
{
  byte data;
  for (int i = 0; i < 8; i ++)
  {
    if (digitalRead (DHpin) == LOW)
    {
      while (digitalRead (DHpin) == LOW); // wait for 50us
      delayMicroseconds (30); // determine the duration of the high level to determine the data is '0 'or '1'
      if (digitalRead (DHpin) == HIGH)
        data |= (1 << (7 - i)); // high front and low in the post
      while (digitalRead (DHpin) == HIGH); // data '1 ', wait for the next one receiver
    }
  }
  return data;
}


void start_test ()
{
  delay (10000);
  digitalWrite (DHpin, LOW); // bus down, send start signal
  delay (30); // delay greater than 18ms, so DHT11 start signal can be detected
  digitalWrite (DHpin, HIGH);
  delayMicroseconds (40); // Wait for DHT11 response delayMicroseconds (40);
  pinMode (DHpin, INPUT);
  while (digitalRead (DHpin) == HIGH);
  delayMicroseconds (80); // DHT11 response, pulled the bus 80us
  while (digitalRead (DHpin) == LOW);
  delayMicroseconds (80); // DHT11 80us after the bus pulled to start delayMicroseconds (80);

  for (int i = 0; i < 4; i++) // receive temperature and humidity for (int i = 0; < 4; ++)
    dat[i] = read_data ();

  pinMode (DHpin, OUTPUT);
  digitalWrite (DHpin, HIGH); // send data once after releasing the bus, digitalWrite (DHpin, HIGH); // send data once after releasing the bus, wait for the host to open next Start signal.
}

void setup ()
{
  delay (10000); 
  
  Serial.begin (9600);
  Serial.flush ();
  pinMode (DHpin, OUTPUT);
}

int z = 9;
void loop ()
{
  if(z==9)
    {
    delay (10000);
    z++;
    }
    
  Serial.print ("Current humdity =");
  Serial.print (dat [0], DEC); // display the humidity Serial.print (dat [0], DEC); //
  Serial.print ('.');
  Serial.print (dat [1], DEC); // display the humidity dec Serial.print (dat [1], DEC); //
  Serial.println ('%');
  Serial.print ("Current temperature =");
  Serial.print (dat [2], DEC); // display the temperature of integer bits;
  Serial.print ('.');
  Serial.print (dat [3], DEC); // display the temperature of decimal
  Serial.println ('C');
  delay (700);

}

