import serial
import csv
import threading
import queue
import time

PORT = '/dev/cu.usbmodem101'
BAUD = 9600
OUTPUT = 'data.csv'

ser = serial.Serial(PORT, BAUD, timeout=1)
ser.reset_input_buffer()
data_queue = queue.Queue()

def reader(ser, q):
    while True:
        try:
            line = ser.readline().decode('utf-8', errors='ignore').strip()
            if line:
                q.put((time.time(), line))
        except Exception:
            break

print(f'Logging to {OUTPUT} ... Ctrl+C to stop')
t = threading.Thread(target=reader, args=(ser, data_queue), daemon=True)
t.start()

count = 0
t0 = None
try:
    with open(OUTPUT, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['time_s', 'pressure_Pa'])
        while True:
            ts, line = data_queue.get()

            # strip optional "something -> " prefix
            idx = line.find(' -> ')
            if idx != -1:
                line = line[idx + 4:]

            # strip unit suffix
            if line.endswith('Pa'):
                line = line[:-2]
            line = line.strip()

            try:
                value = float(line)
            except ValueError:
                continue

            if t0 is None:
                t0 = ts
            writer.writerow([f'{ts - t0:.3f}', value])
            count += 1
            if count % 50 == 0:
                f.flush()
                print(f'{count} samples logged', end='\r')
except KeyboardInterrupt:
    print(f'\nDone. {count} samples saved.')
finally:
    ser.close()